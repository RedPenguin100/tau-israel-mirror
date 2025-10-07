#!/bin/bash
set -e  # Exit on any error
set -x  # Print each command before executing (verbose mode)

echo "=========================================="
echo "=== Render Build Script Starting ==="
echo "=========================================="
echo "Current directory: $(pwd)"
echo "Current user: $(whoami)"
echo "Date: $(date)"
echo ""

# Step 1: Git LFS (skip installation, just pull)
echo "=========================================="
echo "Step 1: Checking Git LFS..."
echo "=========================================="
echo "Git LFS version:"
git lfs version || echo "Git LFS not found, but files may already be pulled during clone"
echo "Ensuring LFS files are pulled..."
git lfs pull || echo "LFS pull skipped (files may already exist)"
echo "✓ Git LFS check complete"
echo "LFS files in repo:"
git lfs ls-files || echo "No LFS tracking or files already pulled"
echo ""

# Step 2: Install Conda (Miniconda)
echo "=========================================="
echo "Step 2: Installing Miniconda..."
echo "=========================================="
echo "Downloading Miniconda installer..."
wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
echo "Installing Miniconda to $HOME/miniconda..."
bash miniconda.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:$PATH"
source $HOME/miniconda/etc/profile.d/conda.sh
echo "Conda version:"
conda --version
echo "✓ Miniconda installed successfully"
echo ""

# Step 3: Create Conda environment from off_target/env.yml
echo "=========================================="
echo "Step 3: Creating Conda environment..."
echo "=========================================="
echo "Checking if off_target/env.yml exists..."
if [ -f "off_target/env.yml" ]; then
    echo "✓ Found off_target/env.yml"
    echo "Contents of env.yml:"
    cat off_target/env.yml
    echo ""
    echo "Creating conda environment 'off_target_env'..."
    conda env create -f off_target/env.yml -n off_target_env
    echo "Activating conda environment..."
    conda activate off_target_env
    echo "Active conda environment:"
    conda info --envs
    echo "Python version in environment:"
    python --version
    echo "✓ Conda environment created and activated"
else
    echo "✗ ERROR: off_target/env.yml not found!"
    exit 1
fi
echo ""

# Step 4: Install pip requirements
echo "=========================================="
echo "Step 4: Installing pip requirements..."
echo "=========================================="
echo "Checking if requirements.txt exists..."
if [ -f "requirements.txt" ]; then
    echo "✓ Found requirements.txt"
    echo "Contents of requirements.txt:"
    cat requirements.txt
    echo ""
    echo "Installing pip packages..."
    pip install -r requirements.txt
    echo "Installed pip packages:"
    pip list
    echo "✓ Pip requirements installed successfully"
else
    echo "✗ ERROR: requirements.txt not found!"
    exit 1
fi
echo ""

echo "=========================================="
echo "=== Build completed successfully! ==="
echo "=========================================="
echo "Build finished at: $(date)"
echo "Total disk usage:"
du -sh $HOME/miniconda 2>/dev/null || echo "Could not calculate disk usage"
echo ""
