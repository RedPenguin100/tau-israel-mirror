#!/bin/bash
set -e

METHOD="h"
MAX_MISMATCHES=3
OUTPUT_PREFIX="res"

# ---------------------------
# Usage / Help
# ---------------------------
usage() {
    cat <<EOF
Usage: $0 [options]
This script performs alignment using Bowtie1 with Hamming distance and creates BAM files.

Options:
  -r <ref.fa>        Reference FASTA file (required)
  -q <reads.fa>      Query FASTA file (required)
  -s <session_id>    Session ID for separating web service sessions (required)
  -o <prefix>        Output prefix (default: res)
  -m <int>           Max mismatches (default: 3)
  -h                 Show this help message and exit

Example:
  $0 -r ref.fa -q reads.fa -s session123 -o results -m 2
EOF
    exit 0
}

# ---------------------------
# Parse command-line arguments
# ---------------------------
while getopts "r:q:s:o:m:h" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        q) QUERY="$OPTARG" ;;
        s) SESSION_ID="$OPTARG" ;;
        o) OUTPUT_PREFIX="$OPTARG" ;;
        m) MAX_MISMATCHES="$OPTARG" ;;
        h) usage ;;
        *) echo "Invalid option"; usage ;;
    esac
done

# ---------------------------
# Validate required arguments
# ---------------------------
for var in REF QUERY SESSION_ID; do
    if [ -z "${!var}" ]; then
        echo "Error: Missing required argument: $var"
        usage
    fi
done

# ---------------------------
# Prepare directory paths
# ---------------------------
REF_BASENAME="${REF##*/}"
REF_BASENAME="${REF_BASENAME%.*}"
QUERY_BASENAME="${QUERY##*/}"
QUERY_BASENAME="${QUERY_BASENAME%.*}"

# Index stays in app directory (reusable across sessions)
INDEX_DIR="index_structure/$REF_BASENAME"

# BAM files go to /tmp (session-specific, temporary)
BAM_DIR="/tmp/bam_files/${SESSION_ID}/${REF_BASENAME}_${QUERY_BASENAME}"

echo "Session ID: $SESSION_ID"
echo "BAM output directory: $BAM_DIR"

# Create the BAM directory in /tmp
mkdir -p "$BAM_DIR"

# ---------------------------
# Build index using build_index.sh
# ---------------------------
echo "Ensuring index is built..."
./build_index.sh -r "$REF" || { echo "Error: build_index.sh failed"; exit 1; }

# ---------------------------
# Alignment function
# ---------------------------
align_bowtie1() {
    SAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sam"
    echo "Running Bowtie1 alignment..."
    bowtie -v "$MAX_MISMATCHES" -a --best -f -x "$INDEX_DIR/${REF_BASENAME}_index" -S "$QUERY" > "$SAM_FILE" \
        || { echo "Error: Bowtie1 alignment failed"; exit 1; }
    echo "SAM file saved: $SAM_FILE"
}

# ---------------------------
# Convert SAM → sorted BAM
# ---------------------------
sam_to_bam() {
    SAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sam"
    BAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sorted.bam"
    echo "Converting SAM to sorted BAM..."
    samtools view -bS "$SAM_FILE" | samtools sort -o "$BAM_FILE" - \
        || { echo "Error: BAM conversion failed"; exit 1; }
    echo "BAM file ready: $BAM_FILE"
}

# ---------------------------
# Index BAM
# ---------------------------
index_bam() {
    BAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sorted.bam"
    if command -v samtools &>/dev/null; then
        echo "Indexing BAM file..."
        samtools index "$BAM_FILE" || { echo "Error: BAM indexing failed"; exit 1; }
        echo "BAM index ready: ${BAM_FILE}.bai"
    else
        echo "Warning: samtools not found. Skipping indexing."
    fi
}

# ---------------------------
# Main workflow
# ---------------------------
if [[ "$METHOD" == "h" || "$METHOD" == "hamming" ]]; then
    align_bowtie1
elif [[ "$METHOD" == "l" || "$METHOD" == "levenstein" ]]; then
    echo "Error: Bowtie2/Levenshtein method not currently supported"
    exit 1
else
    echo "Error: Unknown method: $METHOD"
    exit 1
fi

sam_to_bam
index_bam

echo "Pipeline completed successfully for session: $SESSION_ID"
