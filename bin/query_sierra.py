#!/usr/bin/env python3

# Process all FASTA files

import json
import requests
import os
import sys
import glob
import time
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description='Query Sierra GraphQL API with FASTA sequences')
    parser.add_argument('--query-file', required=True, help='Path to GraphQL query file')
    parser.add_argument('--sierra-port', required=True, help='Port where Sierra server is running')
    parser.add_argument('--fasta-files', nargs='+', required=True, help='FASTA files to process')
    parser.add_argument('--output', default='sierra_results.json', help='Output JSON file')
    parser.add_argument('--individual-dir', default='individual_results', help='Directory for individual results')
    return parser.parse_args()

def read_graphql_query(query_file):
    """Read GraphQL query from file."""
    try:
        with open(query_file, "r") as f:
            graphql_query = f.read().strip()
        print(f"Loaded GraphQL query from {query_file}")
        print(f"Query length: {len(graphql_query)} characters")
        return graphql_query
    except Exception as e:
        print(f"ERROR: Could not read GraphQL query file: {e}", file=sys.stderr)
        sys.exit(1)

def parse_fasta_files(fasta_files):
    """Parse provided FASTA files."""
    print(f"Processing {len(fasta_files)} FASTA files: {fasta_files}")
    
    all_sequences = []
    file_sequence_mapping = {}
    seq_id = 0
    
    for fasta_file in fasta_files:
        print(f"Processing file: {fasta_file}")
        
        if not os.path.exists(fasta_file):
            print(f"ERROR: File not found: {fasta_file}", file=sys.stderr)
            sys.exit(1)
        
        file_sequences = []
        
        try:
            with open(fasta_file, "r") as f:
                content = f.read()
        except Exception as e:
            print(f"ERROR: Could not read {fasta_file}: {e}", file=sys.stderr)
            sys.exit(1)
        
        current_header = None
        current_seq = ""
        
        for line in content.split("\n"):
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequence_entry = {
                        "header": str(seq_id),
                        "sequence": current_seq.replace(" ", "").replace("\n", "")
                    }
                    all_sequences.append(sequence_entry)
                    file_sequences.append(sequence_entry)
                    seq_id += 1
                current_header = line[1:].strip()
                current_seq = ""
            elif line:  # Skip empty lines
                current_seq += line
        
        # Add the last sequence
        if current_header and current_seq:
            sequence_entry = {
                "header": str(seq_id),
                "sequence": current_seq.replace(" ", "").replace("\n", "")
            }
            all_sequences.append(sequence_entry)
            file_sequences.append(sequence_entry)
            seq_id += 1
        
        file_sequence_mapping[fasta_file] = file_sequences
        print(f"  Found {len(file_sequences)} sequences in {fasta_file}")
    
    return all_sequences, file_sequence_mapping


def submit_to_sierra(sierra_url, graphql_query, sequences, max_retries=3):
    """Submit sequences to Sierra GraphQL endpoint with retries."""
    for attempt in range(max_retries):
        try:
            print(f"Submitting sequences to Sierra server (attempt {attempt + 1}/{max_retries})...")
            
            response = requests.post(
                sierra_url,
                json={
                    "query": graphql_query,
                    "variables": {"sequences": sequences}
                },
                headers={"Content-Type": "application/json"},
                timeout=600
            )
            
            print(f"Response status: {response.status_code}")
            
            if response.status_code == 200:
                result_data = response.json()
                
                # Check for GraphQL errors
                if "errors" in result_data:
                    print(f"GraphQL errors: {result_data['errors']}")
                    if attempt == max_retries - 1:
                        sys.exit(1)
                    else:
                        time.sleep(5)
                        continue
                
                return result_data
            else:
                print(f"HTTP Error {response.status_code}: {response.text}")
                if attempt == max_retries - 1:
                    response.raise_for_status()
                else:
                    time.sleep(5)
                    
        except requests.exceptions.RequestException as e:
            print(f"Request error (attempt {attempt + 1}): {e}")
            if hasattr(e, 'response') and e.response is not None:
                print(f"Response content: {e.response.text}")
            
            if attempt == max_retries - 1:
                print(f"ERROR: Failed to submit sequences after {max_retries} attempts", file=sys.stderr)
                sys.exit(1)
            else:
                print(f"Retrying in 10 seconds...")
                time.sleep(10)
        
        except Exception as e:
            print(f"Unexpected error: {e}", file=sys.stderr)
            sys.exit(1)
    
    return None

def save_individual_results(result_data, file_sequence_mapping, output_dir):
    """Save individual results per FASTA file."""
    if "data" not in result_data or "sequenceAnalysis" not in result_data["data"]:
        print("ERROR: Unexpected response structure from Sierra")
        print(f"Response keys: {list(result_data.get('data', {}).keys())}")
        sys.exit(1)
    
    submitted_results = result_data["data"]["sequenceAnalysis"]
    print(f"Received {len(submitted_results)} results from Sierra")
    
    for fasta_file, file_sequences in file_sequence_mapping.items():
        file_results = []
        base_name = os.path.splitext(os.path.basename(fasta_file))[0]
        
        for seq_entry in file_sequences:
            for result in submitted_results:
                if result["inputSequence"]["header"] == seq_entry["header"]:
                    # Create a copy of the result to avoid modifying the original
                    result_copy = result.copy()
                    result_copy["inputSequence"] = result["inputSequence"].copy()
                    result_copy["inputSequence"]["header"] = f"{base_name}"
                    file_results.append(result_copy)
                    break
        
        if file_results:
            individual_file = os.path.join(output_dir, f"{base_name}_sierra.json")
            
            individual_result_data = {
                "data": {
                    "currentVersion": result_data["data"].get("currentVersion"),
                    "currentProgramVersion": result_data["data"].get("currentProgramVersion"),
                    "sequenceAnalysis": file_results
                }
            }
            
            with open(individual_file, "w") as f:
                json.dump(individual_result_data, f, indent=2)
            
            print(f"Saved individual results for {fasta_file} to {individual_file}")
        else:
            print(f"WARNING: No results found for {fasta_file}")

def main():
    args = parse_arguments()
    
    # Read GraphQL query
    graphql_query = read_graphql_query(args.query_file)
    
    # Parse FASTA files
    all_sequences, file_sequence_mapping = parse_fasta_files(args.fasta_files)
    
    print(f"Total sequences to process: {len(all_sequences)}")
    
    if not all_sequences:
        print("ERROR: No sequences found in any FASTA files", file=sys.stderr)
        sys.exit(1)
    
    # Build Sierra URL
    sierra_url = f"http://localhost:{args.sierra_port}/sierra/rest/graphql"
    print(f"Using Sierra URL: {sierra_url}")
    
    # Submit to Sierra
    result_data = submit_to_sierra(sierra_url, graphql_query, all_sequences)
    
    # Save combined results
    with open(args.output, "w") as f:
        json.dump(result_data, f, indent=2)
    
    print(f"Successfully processed {len(all_sequences)} sequences")
    
    # Create directory for individual results
    os.makedirs(args.individual_dir, exist_ok=True)
    
    # Save individual results
    save_individual_results(result_data, file_sequence_mapping, args.individual_dir)



if __name__ == "__main__":
    main()