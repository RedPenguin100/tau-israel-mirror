#!/bin/bash
set -e

METHOD="h"
MAX_MISMATCHES=3
OUTPUT_PREFIX="res"

# Resolve this script's directory so relative helpers work from anywhere
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
BUILD_INDEX="${SCRIPT_DIR}/build_index.sh"

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
  -s <session_id>    Session ID (required)
  -o <prefix>        Output prefix (default: res)
  -m <int>           Max mismatches (default: 3)
  -h                 Show help and exit
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

# TODO: don't hard-code
INDEX_DIR="/tmp/.cache/asodesigner/off_target/index_structure/${REF_BASENAME}"

# BAM files go to /tmp (session-specific)
BAM_DIR="/tmp/bam_files/${SESSION_ID}/${REF_BASENAME}_${QUERY_BASENAME}"

mkdir -p "$BAM_DIR"

# ---------------------------
# Ensure dependencies + index exist
# ---------------------------
command -v bowtie >/dev/null 2>&1 || { echo "Error: bowtie (v1) not found in PATH"; exit 1; }
command -v samtools >/dev/null 2>&1 || { echo "Warning: samtools not found; BAM indexing will be skipped"; }

# build_index.sh via absolute path
[ -f "$BUILD_INDEX" ] || { echo "Error: build_index.sh not found at $BUILD_INDEX"; exit 1; }
[ -x "$BUILD_INDEX" ] || { echo "Error: build_index.sh not executable (chmod +x): $BUILD_INDEX"; exit 1; }

"$BUILD_INDEX" -r "$REF" || { echo "Error: build_index.sh failed"; exit 1; }

# Bowtie expects -x <index_base> (without extensions)
INDEX_BASE="${INDEX_DIR}/${REF_BASENAME}_index"

# ---------------------------
# Alignment function
# ---------------------------
align_bowtie1() {
    SAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sam"
    bowtie -v "$MAX_MISMATCHES" -a --best -f -x "$INDEX_BASE" -S "$QUERY" > "$SAM_FILE" \
        || { echo "Error: Bowtie1 alignment failed"; exit 1; }
}

# ---------------------------
# Convert SAM → sorted BAM
# ---------------------------
sam_to_bam() {
    SAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sam"
    BAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sorted.bam"
    samtools view -bS "$SAM_FILE" | samtools sort -o "$BAM_FILE" - \
        || { echo "Error: BAM conversion failed"; exit 1; }
}

# ---------------------------
# Index BAM
# ---------------------------
index_bam() {
    BAM_FILE="$BAM_DIR/${QUERY_BASENAME}.sorted.bam"
    if command -v samtools &>/dev/null; then
        samtools index "$BAM_FILE" || { echo "Error: BAM indexing failed"; exit 1; }
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
