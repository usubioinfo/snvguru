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

# 5. Environment tool directories
ENV_BIN="$($PKG_MGR run -n "$ENV_NAME" python -c "import sys, os; print(os.path.abspath(os.path.dirname(sys.executable)))")"
ENV_PREFIX="$($PKG_MGR run -n "$ENV_NAME" python -c "import sys, os; print(os.path.abspath(os.path.join(os.path.dirname(sys.executable), '..')))")"

echo -e "${GREEN}[+] Verifying tool binary compatibility...${NC}"

# Check sra-tools compatibility with host GLIBC (e.g. CentOS 7)
if ! $PKG_MGR run -n "$ENV_NAME" prefetch --version &> /dev/null; then
    echo -e "${YELLOW}[!] Host system GLIBC compatibility issue detected for SRA Toolkit (e.g. CentOS 7).${NC}"
    echo -e "${GREEN}[+] Installing NCBI static binaries for SRA Toolkit into environment bin...${NC}"
    TMP_DIR="$(mktemp -d)"
    curl -sL "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-centos_linux64.tar.gz" -o "$TMP_DIR/sratoolkit.tar.gz"
    tar -xzf "$TMP_DIR/sratoolkit.tar.gz" -C "$TMP_DIR"
    cp -rf "$TMP_DIR"/sratoolkit.*-centos_linux64/bin/* "$ENV_BIN/"
    rm -rf "$TMP_DIR"
    echo -e "${GREEN}[+] Static SRA Toolkit configured in $ENV_BIN${NC}"
fi

# Check Trim Galore compatibility
if ! $PKG_MGR run -n "$ENV_NAME" trim_galore --version &> /dev/null; then
    echo -e "${YELLOW}[!] Host system GLIBC compatibility issue detected for trim_galore (e.g. CentOS 7).${NC}"
    echo -e "${GREEN}[+] Installing standalone Trim Galore script into environment bin...${NC}"
    curl -sL "https://raw.githubusercontent.com/FelixKrueger/TrimGalore/0.6.10/trim_galore" -o "$ENV_BIN/trim_galore"
    chmod +x "$ENV_BIN/trim_galore"
    $PKG_MGR run -n "$ENV_NAME" pip install cutadapt
    echo -e "${GREEN}[+] Standalone Trim Galore configured in $ENV_BIN/trim_galore${NC}"
fi

# Check HISAT2 CPU AVX2 instruction compatibility
if ! $PKG_MGR run -n "$ENV_NAME" hisat2-build-s --version &> /dev/null || ! $PKG_MGR run -n "$ENV_NAME" hisat2-align-s --version &> /dev/null; then
    echo -e "${YELLOW}[!] Host CPU compatibility check for hisat2 (e.g. older Xeon without AVX2).${NC}"
    echo -e "${GREEN}[+] Installing generic Linux_x86_64 HISAT2 into environment bin...${NC}"
    $PKG_MGR run -n "$ENV_NAME" python -c "
import urllib.request, zipfile, os, shutil, glob, tempfile
tmp_dir = tempfile.mkdtemp()
zip_path = os.path.join(tmp_dir, 'hisat2.zip')
urllib.request.urlretrieve('https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download', zip_path)
with zipfile.ZipFile(zip_path, 'r') as z:
    z.extractall(tmp_dir)
extracted = glob.glob(os.path.join(tmp_dir, 'hisat2-2.*'))
if extracted:
    for f in glob.glob(os.path.join(extracted[0], 'hisat2*')):
        shutil.copy2(f, '$ENV_BIN')
shutil.rmtree(tmp_dir)
"
    chmod +x "$ENV_BIN"/hisat2*
    echo -e "${GREEN}[+] Generic HISAT2 configured in $ENV_BIN${NC}"
fi

# 6. Ensure JACUSA2 is v2.0.4 (multithreaded LTS release)
JACUSA_JAR_PATH="$ENV_PREFIX/share/jacusa2/JACUSA_v2.0.4.jar"
if [ ! -f "$JACUSA_JAR_PATH" ] || [ ! -f "$ENV_BIN/JACUSA2" ]; then
    echo -e "${GREEN}[+] Setting up JACUSA2 v2.0.4 in environment...${NC}"
    mkdir -p "$(dirname "$JACUSA_JAR_PATH")"
    if [ ! -f "$JACUSA_JAR_PATH" ]; then
        curl -sL "https://github.com/dieterich-lab/JACUSA2/releases/download/v2.0.4/JACUSA_v2.0.4.jar" -o "$JACUSA_JAR_PATH"
    fi
    cat << 'EOF' > "$ENV_BIN/JACUSA2"
#!/usr/bin/env bash
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec "$SCRIPT_DIR/java" -jar "$SCRIPT_DIR/../share/jacusa2/JACUSA_v2.0.4.jar" "$@"
EOF
    chmod +x "$ENV_BIN/JACUSA2"
    echo -e "${GREEN}[+] JACUSA2 v2.0.4 configured in $ENV_BIN/JACUSA2${NC}"
fi

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
