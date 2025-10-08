#!/bin/bash
set -e  # Exit on any error

echo "========================================"
echo "Starting Shell Script Downloads"
echo "========================================"

# Create directory structure
echo "Creating directories..."
mkdir -p /app/aso_gen/data/human/human_v34

# Download genome file
echo ""
echo "[1/1] Downloading GRCh38.p13.genome.fa.gz..."
wget --progress=bar:force:noscroll \
     --timeout=600 \
     --tries=3 \
     https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.p13.genome.fa.gz \
     -O /app/aso_gen/data/human/human_v34/GRCh38.p13.genome.fa.gz

# Unzip the genome file
echo ""
echo "Extracting GRCh38.p13.genome.fa.gz..."
gunzip -f /app/aso_gen/data/human/human_v34/GRCh38.p13.genome.fa.gz

echo ""
echo "✓ Genome file downloaded and extracted!"
echo "✓ Location: /app/aso_gen/data/human/human_v34/GRCh38.p13.genome.fa"
echo "========================================"
