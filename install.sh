#!/usr/bin/env bash

set -e

# ANSI Color codes for clean output
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m' # No Color

echo -e "${BLUE}========================================${NC}"
echo -e "${BLUE}       SNVGuru Automated Installer      ${NC}"
echo -e "${BLUE}========================================${NC}"

# 1. Detect package manager
PKG_MGR=""
if command -v mamba &> /dev/null; then
    PKG_MGR="mamba"
elif command -v micromamba &> /dev/null; then
    PKG_MGR="micromamba"
elif command -v conda &> /dev/null; then
    PKG_MGR="conda"
else
    echo -e "${RED}[ERROR] No Conda package manager detected!${NC}"
    echo -e "Please install Mamba (recommended) or Conda before running this script."
    echo -e "Visit: https://github.com/conda-forge/miniforge#miniforge3"
    exit 1
fi

echo -e "${GREEN}[+] Detected package manager: ${PKG_MGR}${NC}"

# 2. Check if environment 'snvguru' already exists
ENV_NAME="snvguru"
ENV_EXISTS=false

if $PKG_MGR env list | grep -E -q "^[[:space:]]*${ENV_NAME}([[:space:]]|$)"; then
    ENV_EXISTS=true
fi

# 3. Create or update the Conda environment
if [ "$ENV_EXISTS" = true ]; then
    echo -e "${YELLOW}[!] Environment '${ENV_NAME}' already exists. Updating environment...${NC}"
    $PKG_MGR env update -n "$ENV_NAME" -f environment.yml --prune
else
    echo -e "${GREEN}[+] Creating environment '${ENV_NAME}' from environment.yml...${NC}"
    $PKG_MGR env create -f environment.yml
fi

# 4. Install the package in editable mode inside the environment
echo -e "${GREEN}[+] Installing SNVGuru package into '${ENV_NAME}'...${NC}"
$PKG_MGR run -n "$ENV_NAME" pip install -e .

echo ""
echo -e "${GREEN}====================================================${NC}"
echo -e "${GREEN}      Installation completed successfully!          ${NC}"
echo -e "${GREEN}====================================================${NC}"
echo ""
echo -e "To start using SNVGuru:"
echo -e "  1. Activate the environment:"
echo -e "     ${BLUE}${PKG_MGR} activate ${ENV_NAME}${NC}"
echo ""
echo -e "  2. Go to your working directory and initialize your project:"
echo -e "     ${BLUE}snvguru init${NC}"
echo ""
echo -e "  3. Run the pipeline:"
echo -e "     ${BLUE}snvguru${NC}"
echo ""
