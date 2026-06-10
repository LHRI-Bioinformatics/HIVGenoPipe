#!/usr/bin/env python3

import sys
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from collections import Counter

# --- Default deletion-window parameters (also passable as CLI args 4 & 5) ---
DEFAULT_MIN_DEL = 0   # shorter runs -> re-insert reference (frameshift/codon guard)
DEFAULT_MAX_DEL = 6   # longer runs  -> re-insert reference (reads portray it anyway)


def read_sequences_from_fasta(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(str(record.seq))
    return sequences


def write_sequence_to_fasta(id, sequence, file_path, description=""):
    seq_record = SeqRecord(Seq(sequence), id=id, description=description)
    with open(file_path, 'w') as output_handle:
        SeqIO.write(seq_record, output_handle, 'fasta')


def contig_coverage_bounds(seq):
    """Return (first, last) column indices a contig actually covers.

    Leading/trailing gaps are NOT counted as coverage — the contig simply
    hasn't reached that region and should not cast a phantom deletion vote.
    Internal gaps within [first, last] DO count and are treated as real
    deletion votes. Returns None for an all-gap (empty) sequence.
    """
    stripped = seq.strip("-")
    if not stripped:
        return None
    first = len(seq) - len(seq.lstrip("-"))
    last  = len(seq.rstrip("-")) - 1
    return first, last


def column_call(column_chars, ref_base):
    """Plurality call for a column, INCLUDING '-' as a valid vote.

    Args:
        column_chars : characters from contigs that COVER this column
                       (may include '-' for internal deletions).
        ref_base     : reference character, used only as a tiebreaker.

    Returns one of:
        a base / '-'  — clear plurality winner
        ref_base      — tie broken by reference (ref is one of the tied chars)
        sorted pick   — tie with no reference match; deterministic fallback
        None          — no contig covers this column at all
    """
    if not column_chars:
        return None

    counts = Counter(column_chars)
    top_count = counts.most_common(1)[0][1]
    top_chars = [c for c, n in counts.items() if n == top_count]

    if len(top_chars) == 1:
        return top_chars[0]
    elif ref_base in top_chars:
        return ref_base          # reference breaks the tie
    else:
        return sorted(top_chars)[0]   # deterministic fallback (no ref match)


def create_hybrid_consensus(sequences, min_del=DEFAULT_MIN_DEL, max_del=DEFAULT_MAX_DEL):
    """Build a hybrid consensus from a MAFFT/PAGAN MSA.

    sequences[0] is the HXB2 reference; sequences[1:] are assembled contigs.

    Decision rules applied at each alignment column (in priority order):
      1. HONORED DELETION — output '-' (stripped at end) when ALL of:
            a. the plurality of covering contigs voted '-'
            b. the run of consecutive deletion columns is min_del..max_del long
            c. the run is internal (contig coverage exists on both flanks)
         Outside the length window, or at edges, the reference is re-inserted.
      2. CONTIG MAJORITY — among contigs that cover this column, take the
         plurality base (gaps now count as votes). Reference breaks ties;
         if ref is not a tied candidate, take sorted(tied_chars)[0].
      3. NO COVERAGE FALLBACK — if no contig covers this column, use the
         reference base (preserves HXB2 coordinates at ends and in large gaps).
    """
    if not all(len(s) == len(sequences[0]) for s in sequences):
        raise ValueError("All sequences must be of the same length")

    reference = str(sequences[0])
    aln_len   = len(reference)
    contigs   = sequences[1:]

    # ---- Reference flank trimming bookkeeping (unchanged from original) ----
    number_right_strip_bases = len(reference) - len(reference.rstrip("-"))
    number_left_strip_bases  = len(reference) - len(reference.lstrip("-"))
    print("right stripped bases =", number_right_strip_bases)
    print("left strip bases     =", number_left_strip_bases)

    # ---- Per-contig covered spans (leading/trailing gaps excluded) ----
    bounds = [contig_coverage_bounds(c) for c in contigs]

    # =====================================================================
    # PASS 1 — per-column consensus call with '-' as a valid vote.
    #   resolved[i]    : winning character (base or '-') at column i
    #   has_coverage[i]: True if at least one contig covers column i
    # =====================================================================
    resolved     = []
    has_coverage = []

    for i in range(aln_len):
        column_chars = []
        for c, b in zip(contigs, bounds):
            if b is None:
                continue
            first, last = b
            if first <= i <= last:
                column_chars.append(c[i])   # base or internal '-'

        call = column_call(column_chars, reference[i])

        if call is None:
            resolved.append(reference[i])   # Rule 3: no coverage -> reference
            has_coverage.append(False)
        else:
            resolved.append(call)
            has_coverage.append(True)

    # =====================================================================
    # PASS 2 — walk maximal '-' runs and decide: HONOR or REVERT.
    #   Honor  => keep '-' (stripped at assembly)
    #   Revert => paste reference bases back across the run
    # =====================================================================
    consensus_chars = list(resolved)

    i = 0
    while i < aln_len:
        if consensus_chars[i] != '-':
            i += 1
            continue

        run_start = i
        while i < aln_len and consensus_chars[i] == '-':
            i += 1
        run_end = i   # exclusive
        run_len = run_end - run_start

        contig_supported = any(has_coverage[j] for j in range(run_start, run_end))
        covered_before   = any(has_coverage[j] for j in range(0, run_start))
        covered_after    = any(has_coverage[j] for j in range(run_end, aln_len))
        is_internal      = covered_before and covered_after

        honor = (
            contig_supported
            and is_internal
            and min_del <= run_len <= max_del
        )

        if honor:
            print("HONORED deletion  : cols {:>6}-{:>6}  len={}".format(
                run_start, run_end - 1, run_len))
        else:
            reason = (
                "no-contig-support" if not contig_supported else
                "flank/edge"        if not is_internal      else
                "len<{}".format(min_del) if run_len < min_del else
                "len>{}".format(max_del)
            )
            print("REVERTED to ref   : cols {:>6}-{:>6}  len={}  ({})".format(
                run_start, run_end - 1, run_len, reason))
            for j in range(run_start, run_end):
                consensus_chars[j] = reference[j]

    # =====================================================================
    # Final assembly: drop honored '-' columns, then trim reference flanks.
    # =====================================================================
    hybrid_consensus = ''.join(ch for ch in consensus_chars if ch != '-')

    if number_left_strip_bases != 0:
        hybrid_consensus = hybrid_consensus[number_left_strip_bases:]
    if number_right_strip_bases != 0:
        hybrid_consensus = hybrid_consensus[:-number_right_strip_bases]

    print(">final hybrid\n", hybrid_consensus)
    print("length of hybrid =", len(hybrid_consensus))

    return hybrid_consensus


def main():
    alignment_file    = sys.argv[1]   # MAFFT/PAGAN MSA FASTA
    hybrid_header     = sys.argv[2]   # output sequence ID
    hybrid_output_file = sys.argv[3]  # output FASTA path

    # Optional overrides for the deletion-length window
    min_del = int(sys.argv[4]) if len(sys.argv) > 4 else DEFAULT_MIN_DEL
    max_del = int(sys.argv[5]) if len(sys.argv) > 5 else DEFAULT_MAX_DEL

    aligned_sequences = read_sequences_from_fasta(alignment_file)
    hybrid_consensus  = create_hybrid_consensus(aligned_sequences, min_del=min_del, max_del=max_del)
    write_sequence_to_fasta(hybrid_header, hybrid_consensus, hybrid_output_file)


if __name__ == "__main__":
    main()