#!/bin/bash
# aso_counts.sh - Modified for better error reporting to Python subprocess

# Resolve this script's directory for reliable relative paths
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
GET_ASO="${SCRIPT_DIR}/get_aso_off_target.sh"
CAT_PY="${SCRIPT_DIR}/cat.py"
COUNTER_PY="${SCRIPT_DIR}/counter.py"

# Send errors to stdout so Python can capture them
error_exit() {
    echo "__ERROR__: $1"
    echo "__ERROR_LINE__: $2"
    echo "__ERROR_CMD__: $3"
    exit 1
}

trap 'error_exit "Command failed" $LINENO "$BASH_COMMAND"' ERR
set -e

# -------------------------------------------
# Parse arguments
# -------------------------------------------
while getopts "q:k:s:t:h" opt; do
  case $opt in
    q) ASO_QUERY=$OPTARG ;;
    k) K=$OPTARG ;;
    s) SESSION_ID=$OPTARG ;;
    t) TARGET_FILE=$OPTARG ;;
    h)
      echo "Usage: $0 -q aso_query_sequence -k parameter_K_value -s session_id [-t target_file]"
      exit 0
      ;;
    *)
      error_exit "Invalid flag" $LINENO "flag parsing"
      ;;
  esac
done

# Validate required arguments
if [[ -z "$ASO_QUERY" ]]; then
  error_exit "Missing -q (ASO_QUERY)" $LINENO "arg check"
fi

if [[ -z "$K" ]]; then
  error_exit "Missing -k (K value)" $LINENO "arg check"
fi

if [[ -z "$SESSION_ID" ]]; then
  error_exit "Missing -s (SESSION_ID)" $LINENO "arg check"
fi

# Validate files exist
if [[ ! -f "$ASO_QUERY" ]]; then
  error_exit "ASO query file not found: $ASO_QUERY" $LINENO "file check"
fi

if [[ -n "$TARGET_FILE" && ! -f "$TARGET_FILE" ]]; then
  error_exit "Target file not found: $TARGET_FILE" $LINENO "file check"
fi

# Validate helper scripts exist
[[ ! -f "$GET_ASO" ]] && error_exit "Helper not found: $GET_ASO" $LINENO "file check"
[[ ! -x "$GET_ASO" ]] && error_exit "Helper not executable: $GET_ASO (chmod +x)" $LINENO "perm check"
[[ ! -f "$CAT_PY" ]] && error_exit "cat.py not found: $CAT_PY" $LINENO "file check"
[[ ! -f "$COUNTER_PY" ]] && error_exit "counter.py not found: $COUNTER_PY" $LINENO "file check"

#echo "__INFO__: Starting pipeline for session $SESSION_ID"
#echo "__INFO__: ASO_QUERY=$ASO_QUERY, K=$K, TARGET_FILE=${TARGET_FILE:-none}"
#echo "__INFO__: SCRIPT_DIR=$SCRIPT_DIR"

# Run pipeline
if [[ -n "$TARGET_FILE" ]]; then
#  echo "__INFO__: Using target file: $TARGET_FILE"
  "$GET_ASO" -r "$TARGET_FILE" -q "$ASO_QUERY" -k "$K" -s "$SESSION_ID"
else
  chromosomes="/tmp/.cache/asodesigner/chromosomes"

  if [[ ! -d "$chromosomes" ]]; then
    error_exit "Chromosomes directory not found: $chromosomes" $LINENO "dir check"
  fi

  chr_count=$(ls -1 "$chromosomes"/*.fa 2>/dev/null | wc -l)
  if [[ $chr_count -eq 0 ]]; then
    error_exit "No .fa files in $chromosomes" $LINENO "file count"
  fi

#  echo "__INFO__: Processing $chr_count chromosome files"

  for fa_file in "$chromosomes"/*.fa; do
    [[ ! -f "$fa_file" ]] && continue
#    echo "__INFO__: Processing $fa_file"
    "$GET_ASO" -r "$fa_file" -q "$ASO_QUERY" -k "$K" -s "$SESSION_ID"
  done
fi

# Extract ASO basename
ASO_BASENAME=$(basename "$ASO_QUERY")
ASO_BASENAME=${ASO_BASENAME%.fa}
ASO_BASENAME=${ASO_BASENAME%.fasta}

# Work in session results directory
SESSION_RES_DIR="/tmp/res/${SESSION_ID}"
mkdir -p "$SESSION_RES_DIR"
cd "$SESSION_RES_DIR"

#echo "__INFO__: Running cat.py"
python "$CAT_PY" -d .

if [[ ! -f "genome_hits.json" ]]; then
  error_exit "cat.py did not create genome_hits.json" $LINENO "cat.py output"
fi

echo "__INFO__: Running counter.py"
python "$COUNTER_PY"

if [[ ! -f "genome_hits_count.json" ]]; then
  error_exit "counter.py did not create genome_hits_count.json" $LINENO "counter.py output"
fi

# Rename output
OUTPUT_FILE="genome_hits_count_${ASO_BASENAME}.json"
mv "genome_hits_count.json" "$OUTPUT_FILE"

if [[ ! -f "$OUTPUT_FILE" ]]; then
  error_exit "Output file missing after rename: $OUTPUT_FILE" $LINENO "final check"
fi

#echo "__INFO__: Pipeline completed successfully"
echo "__RETURN__: $(pwd)/$OUTPUT_FILE"
