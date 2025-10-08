#!/bin/bash
set -e

# ---------------------------
# Default parameters
# ---------------------------
TEST_MODE=0

# ---------------------------
# Usage / Help
# ---------------------------
usage() {
    echo "Usage: $0 -r <reference.fa> -q <query.fa> -k <max_hamming_distance> -s <session_id> [-t]"
    echo
    echo "  -r    Reference file (FASTA)"
    echo "  -q    Query file (FASTA)"
    echo "  -k    Maximum Hamming distance allowed"
    echo "  -s    Session ID for separating web service sessions (required)"
    echo "  -t    Test mode (does not execute main commands)"
    echo "  -h    Show this help message and exit"
    exit 1
}

# ---------------------------
# Parse arguments
# ---------------------------
while getopts "r:q:k:s:th" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        q) QUERY="$OPTARG" ;;
        k) K="$OPTARG" ;;
        s) SESSION_ID="$OPTARG" ;;
        t) TEST_MODE=1 ;;
        h) usage ;;
        *) usage ;;
    esac
done

# ---------------------------
# Check required args
# ---------------------------
if [[ $TEST_MODE -eq 0 ]]; then
    if [[ -z "$REF" || -z "$QUERY" || -z "$K" || -z "$SESSION_ID" ]]; then
        echo "Error: Missing required arguments."
        usage
    fi
else
    REF="fasta_files/ref.fasta"
    QUERY="fasta_files/aso_reads.1.fasta"
    K=2
    SESSION_ID="test_session"
fi

# ---------------------------
# Determine output prefix from reference and query
# ---------------------------
REF_BASENAME="${REF##*/}"
REF_BASENAME="${REF_BASENAME%.*}"
QUERY_BASENAME="${QUERY##*/}"
QUERY_BASENAME="${QUERY_BASENAME%.*}"

# Auto-generate output subdirectory with session ID in /tmp
FIND_OUTDIR="${REF_BASENAME}_${QUERY_BASENAME}"
BAM_DIR="/tmp/bam_files/${SESSION_ID}/${REF_BASENAME}_${QUERY_BASENAME}"
BAM_FILE="${BAM_DIR}/${QUERY_BASENAME}.sorted.bam"

echo "Session ID: $SESSION_ID"
echo "BAM file: $BAM_FILE"

# Output directory for find_matches_in_bam.py with session ID in /tmp
RES_DIR="/tmp/res/${SESSION_ID}/${FIND_OUTDIR}"

# ---------------------------
# Run main commands
# ---------------------------
./create_bam.sh -r "$REF" -q "$QUERY" -s "$SESSION_ID" -o "$REF_BASENAME" -m "$K"

if [[ $TEST_MODE -eq 1 ]]; then
    mkdir -p /tmp/test
    ./find_matches_in_bam.py -b "$BAM_FILE" -s - -mm "$K" -o /tmp/test 
else
    mkdir -p "$RES_DIR"
    ./find_matches_in_bam.py -b "$BAM_FILE" -s - -mm "$K" -o "$RES_DIR" --single-json
fi

echo "Pipeline completed successfully for session: $SESSION_ID"
