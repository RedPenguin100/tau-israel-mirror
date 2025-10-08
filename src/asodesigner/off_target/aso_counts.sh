#!/bin/bash

# -------------------------------------------
# Helper (-h)
# -------------------------------------------
# Usage:
#   ./run_get_aso_off_target.sh -q ASO_QUERY -k K -s SESSION_ID [-t TARGET_FILE]
#
# Arguments:
#   -q  aso_query_sequence (FASTA file)
#   -k  parameter_K_value
#   -s  session_id (required for web service)
#   -t  target_file (optional, if not provided uses chromosomes directory)
# -------------------------------------------
while getopts "q:k:s:t:h" opt; do
  case $opt in
    q) ASO_QUERY=$OPTARG ;;
    k) K=$OPTARG ;;
    s) SESSION_ID=$OPTARG ;;
    t) TARGET_FILE=$OPTARG ;;
    h)
      echo "Usage: $0 -q aso_query_sequence -k parameter_K_value -s session_id [-t target_file]"
      echo ""
      echo "If -t is provided, uses that target file."
      echo "If -t is not provided, processes all files in ./chromosomes/ directory."
      exit 0
      ;;
    *)
      echo "Invalid flag. Use -h for help."
      exit 1
      ;;
  esac
done

# Check required args
if [[ -z "$ASO_QUERY" || -z "$K" || -z "$SESSION_ID" ]]; then
  echo "Error: Missing required arguments. Use -h for help."
  exit 1
fi

echo "Starting pipeline for session: $SESSION_ID"

# Check if target file was provided
if [[ -n "$TARGET_FILE" ]]; then
  # Use the provided target file
  echo "Using target file: $TARGET_FILE"
  ./get_aso_off_target.sh -r "$TARGET_FILE" -q "$ASO_QUERY" -k "$K" -s "$SESSION_ID"
else
  # Use all chromosome files
  chromosomes="../chromosomes"
  echo "Processing all chromosome files from: $chromosomes"
  
  for fa_file in "$chromosomes"/*.fa; do
    echo "Processing $fa_file ..."
    ./get_aso_off_target.sh -r "$fa_file" -q "$ASO_QUERY" -k "$K" -s "$SESSION_ID"
  done
fi

# Extract ASO basename
ASO_BASENAME=$(basename "$ASO_QUERY")
ASO_BASENAME=${ASO_BASENAME%.fa}
ASO_BASENAME=${ASO_BASENAME%.fasta}

# Work in session-specific results directory in /tmp
SESSION_RES_DIR="/tmp/res/${SESSION_ID}"
CODE_DIR=$(pwd)
mkdir -p "$SESSION_RES_DIR"
cd "$SESSION_RES_DIR"

# Run Python processing scripts
python "${CODE_DIR}/cat.py" -d .
python "${CODE_DIR}/counter.py"

# Rename output with ASO basename
OUTPUT_FILE="genome_hits_count_${ASO_BASENAME}.json"
mv "genome_hits_count.json" "$OUTPUT_FILE"

echo "__RETURN__: $(pwd)/$OUTPUT_FILE"
echo "Pipeline completed successfully for session: $SESSION_ID"
