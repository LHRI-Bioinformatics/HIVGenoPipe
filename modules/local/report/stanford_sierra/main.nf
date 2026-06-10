process REPORT_STANFORD_SIERRA {
    tag "$meta.id"
    
    input:
    tuple val(meta), path(fastas)
    path graphql_query
    
    output:
    // path "sierra_results.json", emit: json
    tuple val(meta), path ("individual_results/*.json"), emit: individual_json, optional: true
    tuple val(meta), path ("individual_results/*Amb*.json"), emit: amb_jsons, optional: true
    tuple val(meta), path ("individual_results/*consensus*.json"), emit: consensus_json, optional: true
    path("sierra_results.json"),  emit: combined_json, optional: true
    path "versions.yml", emit: versions
    
    script:
    """
    # Find an available port
    SIERRA_PORT=\$(python3 -c "import socket; s=socket.socket(); s.bind(('',0)); print(s.getsockname()[1]); s.close()")
    echo "Using port: \$SIERRA_PORT"
    
    # Generate unique container name
    CONTAINER_NAME="sierra-\${RANDOM}-\$\$"
    echo "Container name: \$CONTAINER_NAME"
    
    # Start Sierra container with correct port mapping
    # Container runs on 8080 internally, we map it to our dynamic port externally
    echo "Starting Sierra container..."
    echo "Port mapping: \$SIERRA_PORT:8080 (host:container)"
    docker run -d --name \$CONTAINER_NAME -p \$SIERRA_PORT:8080 hivdb/sierra:latest dev
    
    # Wait for container to be created
    sleep 3
    
    CONTAINER_ID=\$(docker ps -q --filter "name=\$CONTAINER_NAME" | head -1)
    
    if [ -z "\$CONTAINER_ID" ]; then
        echo "ERROR: Failed to start Sierra container" >&2
        docker logs \$CONTAINER_NAME 2>&1 || true
        exit 1
    fi
    
    echo "Started Sierra container: \$CONTAINER_ID"
    echo "Sierra GraphQL endpoint: http://localhost:\$SIERRA_PORT/sierra/rest/graphql"
    
    # Wait for Sierra to be fully ready
    echo "Waiting for Sierra server to be ready (this may take 2-3 minutes)..."
    
    for i in {1..180}; do  # 3 minutes timeout
        # Check container is still running
        if ! docker ps -q --filter "id=\$CONTAINER_ID" | grep -q \$CONTAINER_ID; then
            echo "ERROR: Container stopped unexpectedly"
            docker logs \$CONTAINER_ID 2>&1
            exit 1
        fi
        
        # Check if we can connect to the port and GraphQL endpoint is ready
        if nc -z localhost \$SIERRA_PORT 2>/dev/null || curl -s --connect-timeout 1 http://localhost:\$SIERRA_PORT/ >/dev/null 2>&1; then
            # Try the Sierra GraphQL endpoint
            GRAPHQL_STATUS=\$(curl -s -o /dev/null -w "%{http_code}" http://localhost:\$SIERRA_PORT/sierra/rest/graphql 2>/dev/null || echo "000")
            
            # Accept 200, 400, or 405 as "ready"
            if [ "\$GRAPHQL_STATUS" = "200" ] || [ "\$GRAPHQL_STATUS" = "400" ] || [ "\$GRAPHQL_STATUS" = "405" ]; then
                echo "Sierra server appears ready!"
                break
            fi
        fi
        
        if [ \$i -eq 180 ]; then
            echo "ERROR: Sierra server failed to start within 3 minutes"
            echo "Final container logs:"
            docker logs \$CONTAINER_ID 2>&1 | tail -50
            docker stop \$CONTAINER_ID 2>/dev/null || true
            docker rm \$CONTAINER_ID 2>/dev/null || true
            exit 1
        fi
        
        # Show progress every 15 seconds
        if [ \$((i % 15)) -eq 0 ]; then
            echo "Still waiting... (\$i/180 seconds)"
        fi
        
        sleep 1
    done
    
    echo "Sierra server is ready on port \$SIERRA_PORT!"
    
    # Test the GraphQL endpoint
    echo "Testing GraphQL endpoint..."
    TEST_RESPONSE=\$(curl -s -X POST \\
        -H "Content-Type: application/json" \\
        -d '{"query": "query { viewer { currentVersion { text } } }"}' \\
        http://localhost:\$SIERRA_PORT/sierra/rest/graphql 2>&1)
    
    echo "Test response: \$TEST_RESPONSE"
    
    
    # Create directory for individual results
    mkdir -p individual_results
    
    query_sierra.py \\
        --query-file ${graphql_query} \\
        --sierra-port \$SIERRA_PORT \\
        --fasta-files ${fastas} \\
        --output sierra_results.json \\
        --individual-dir individual_results
    
    # Clean up container
    echo "Cleaning up Sierra container..."
    docker stop \$CONTAINER_ID 2>/dev/null || true
    docker rm \$CONTAINER_ID 2>/dev/null || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sierra: "latest"
        docker: \$(docker --version | cut -d' ' -f3 | tr -d ',')
    END_VERSIONS
    """
}