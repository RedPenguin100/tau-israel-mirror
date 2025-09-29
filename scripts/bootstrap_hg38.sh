#!/usr/bin/env bash
set -euo pipefail

METHOD=${GENOME_FETCH_METHOD:-curl}
TARGET_DIR=${HG38_TARGET_DIR:-/genome-cache}
URL=${ASODESIGNER_HG38_URL:-https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz}
ARCHIVE_PATH="${TARGET_DIR}/hg38.fa.gz"
FASTA_PATH="${TARGET_DIR}/hg38.fa"

mkdir -p "${TARGET_DIR}"

if [ -s "${FASTA_PATH}" ]; then
  echo "hg38 already present at ${FASTA_PATH}"
  exit 0
fi

case "${METHOD}" in
  curl)
    echo "Downloading hg38 with curl"
    curl -L "${URL}" -o "${ARCHIVE_PATH}"
    ;;
  python)
    echo "Downloading hg38 with Python requests"
    python3 /workspace/scripts/download_hg38_requests.py "${URL}" "${ARCHIVE_PATH}"
    ;;
  *)
    echo "Unknown GENOME_FETCH_METHOD: ${METHOD}" >&2
    exit 1
    ;;
esac

if [ ! -s "${ARCHIVE_PATH}" ]; then
  echo "Download failed" >&2
  exit 1
fi

gunzip -c "${ARCHIVE_PATH}" > "${FASTA_PATH}"
echo "hg38 ready at ${FASTA_PATH}"
