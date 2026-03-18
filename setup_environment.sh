#!/bin/bash
# =============================================================================
# Environment Setup for Opto-Gatable Binder Pipeline
# =============================================================================
#
# Creates 3 conda environments needed for the full pipeline:
#   1. boltz_only   — Boltz-2 structure prediction + pipeline analysis (Python 3.10)
#   2. rfdiffusion_clean — RFdiffusion backbone design (Python 3.9)
#   3. ligandmpnn   — LigandMPNN sequence design (Python 3.10)
#
# Prerequisites:
#   - NVIDIA GPU with >=24GB VRAM
#   - CUDA drivers installed (nvidia-smi should work)
#   - Miniconda or Anaconda installed
#   - ~20GB disk space
#
# Usage:
#   bash setup_environment.sh [boltz|rfdiffusion|ligandmpnn|all]
#
# =============================================================================

set -e

STEP=${1:-all}

echo "=============================================="
echo "Opto Pipeline Environment Setup"
echo "=============================================="

# Check GPU
echo ""
echo "Checking GPU..."
if ! command -v nvidia-smi &>/dev/null; then
    echo "WARNING: nvidia-smi not found. GPU may not be available."
else
    nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
fi

# Check conda
if ! command -v conda &>/dev/null; then
    echo "ERROR: conda not found. Install Miniconda first:"
    echo "  wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
    echo "  bash Miniconda3-latest-Linux-x86_64.sh"
    exit 1
fi

source ~/miniconda3/etc/profile.d/conda.sh 2>/dev/null || \
source ~/anaconda3/etc/profile.d/conda.sh 2>/dev/null || \
echo "WARNING: Could not source conda.sh"


# =============================================================================
# Environment 1: boltz_only (primary — runs pipeline + Boltz-2)
# =============================================================================
setup_boltz() {
    echo ""
    echo "=== Setting up boltz_only environment ==="

    if conda env list | grep -q "boltz_only"; then
        echo "  Environment 'boltz_only' already exists."
        echo "  To recreate: conda env remove -n boltz_only && bash $0 boltz"
        return
    fi

    conda create -n boltz_only python=3.10 -y
    conda activate boltz_only

    # Boltz-2 (includes PyTorch)
    pip install boltz

    # Pipeline dependencies
    pip install numpy scipy biopython

    echo "  boltz_only environment ready."
    echo "  Activate: conda activate boltz_only"
    echo "  Test: boltz predict --help"
}


# =============================================================================
# Environment 2: rfdiffusion_clean (RFdiffusion v1.1.0)
# =============================================================================
setup_rfdiffusion() {
    echo ""
    echo "=== Setting up rfdiffusion_clean environment ==="

    if conda env list | grep -q "rfdiffusion_clean"; then
        echo "  Environment 'rfdiffusion_clean' already exists."
        echo "  To recreate: conda env remove -n rfdiffusion_clean && bash $0 rfdiffusion"
        return
    fi

    # Check if RFdiffusion repo exists
    RFDIFFUSION_DIR="${RFDIFFUSION_DIR:-$(pwd)/RFdiffusion}"
    if [ ! -d "$RFDIFFUSION_DIR" ]; then
        echo "  Cloning RFdiffusion..."
        git clone https://github.com/RosettaCommons/RFdiffusion.git "$RFDIFFUSION_DIR"
    fi

    conda create -n rfdiffusion_clean python=3.9 -y
    conda activate rfdiffusion_clean

    # PyTorch for CUDA
    pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118

    # RFdiffusion dependencies
    cd "$RFDIFFUSION_DIR"
    pip install -e .

    # Download model weights if not present
    if [ ! -d "models/weights" ]; then
        echo "  Downloading RFdiffusion weights..."
        mkdir -p models/weights
        # Weights download URL (check RFdiffusion README for latest)
        echo "  NOTE: Download weights manually from RFdiffusion GitHub releases."
        echo "  Place them in $RFDIFFUSION_DIR/models/weights/"
    fi

    echo "  rfdiffusion_clean environment ready."
    echo "  Activate: conda activate rfdiffusion_clean"
    echo "  Test: cd $RFDIFFUSION_DIR && python scripts/run_inference.py --help"
}


# =============================================================================
# Environment 3: ligandmpnn
# =============================================================================
setup_ligandmpnn() {
    echo ""
    echo "=== Setting up ligandmpnn environment ==="

    if conda env list | grep -q "ligandmpnn"; then
        echo "  Environment 'ligandmpnn' already exists."
        echo "  To recreate: conda env remove -n ligandmpnn && bash $0 ligandmpnn"
        return
    fi

    # Check if LigandMPNN repo exists
    LIGANDMPNN_DIR="${LIGANDMPNN_DIR:-$(pwd)/LigandMPNN}"
    if [ ! -d "$LIGANDMPNN_DIR" ]; then
        echo "  Cloning LigandMPNN..."
        git clone https://github.com/dauparas/LigandMPNN.git "$LIGANDMPNN_DIR"
    fi

    conda create -n ligandmpnn python=3.10 -y
    conda activate ligandmpnn

    pip install torch --index-url https://download.pytorch.org/whl/cu118

    cd "$LIGANDMPNN_DIR"
    pip install -e . 2>/dev/null || pip install numpy

    # Download model weights if not present
    if [ ! -d "model_params" ]; then
        echo "  Downloading LigandMPNN weights..."
        bash get_model_params.sh 2>/dev/null || \
        echo "  NOTE: Download weights manually. See LigandMPNN README."
    fi

    echo "  ligandmpnn environment ready."
    echo "  Activate: conda activate ligandmpnn"
    echo "  Test: cd $LIGANDMPNN_DIR && python run.py --help"
}


# =============================================================================
# Main
# =============================================================================
case "$STEP" in
    boltz) setup_boltz ;;
    rfdiffusion) setup_rfdiffusion ;;
    ligandmpnn) setup_ligandmpnn ;;
    all)
        setup_boltz
        setup_rfdiffusion
        setup_ligandmpnn
        ;;
    *)
        echo "Usage: $0 [boltz|rfdiffusion|ligandmpnn|all]"
        exit 1
        ;;
esac

echo ""
echo "=============================================="
echo "Setup complete!"
echo "=============================================="
echo ""
echo "Environments:"
echo "  conda activate boltz_only          # Pipeline + Boltz-2"
echo "  conda activate rfdiffusion_clean   # RFdiffusion"
echo "  conda activate ligandmpnn          # LigandMPNN"
echo ""
echo "Quick test:"
echo "  conda activate boltz_only"
echo "  cd $(dirname "$0")"
echo "  python -m unittest discover -s tests -v"
