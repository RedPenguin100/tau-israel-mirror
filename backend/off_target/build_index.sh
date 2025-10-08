#!/bin/bash
set -e

# ---------------------------
# Usage / Help
# ---------------------------
usage() {
    cat <<EOF
Usage: $0 [options]
This script builds a Bowtie1 index for a reference FASTA file.

Options:
  -r <ref.fa>        Reference FASTA file (required)
  -h                 Show this help message and exit

Example:
  $0 -r ref.fa
EOF
    exit 0
}

# ---------------------------
# Parse command-line arguments
# ---------------------------
while getopts "r:h" opt; do
    case $opt in
        r) REF="$OPTARG" ;;
        h) usage ;;
        *) echo "Invalid option"; usage ;;
    esac
done

# ---------------------------
# Validate required arguments
# ---------------------------
if [ -z "$REF" ]; then
    echo "Error: Missing required argument -r <ref.fa>"
    usage
fi

# ---------------------------
# Prepare directory paths
# ---------------------------
REF_BASENAME="${REF##*/}"           # remove path
REF_BASENAME="${REF_BASENAME%.*}"   # remove extension after the last dot
INDEX_DIR="index_structure/$REF_BASENAME"

# ---------------------------
# Build Bowtie1 index
# ---------------------------
if [ -d "$INDEX_DIR" ]; then
    echo "Index directory exists: $INDEX_DIR"
    echo "Skipping index building."
else
    mkdir -p "$INDEX_DIR"
    echo "Building Bowtie1 index for: $REF"
    bowtie-build "$REF" "$INDEX_DIR/${REF_BASENAME}_index" \
        || { echo "Error: bowtie-build failed"; exit 1; }
    echo "Index built successfully: $INDEX_DIR/${REF_BASENAME}_index"
fi
