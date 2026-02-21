#!/bin/bash
# Download and compile PHYLIP 3.697 for validation testing.
#
# This script downloads the official PHYLIP distribution from
# Joe Felsenstein's website and compiles it on macOS/Linux.
#
# After running this script, the PHYLIP executables will be available
# at validation/phylip-3.697/exe/
#
# Usage:
#   cd validation && bash setup.sh

set -euo pipefail

PHYLIP_VERSION="3.697"
PHYLIP_DIR="phylip-${PHYLIP_VERSION}"
PHYLIP_TARBALL="${PHYLIP_DIR}.tar.gz"
PHYLIP_URL="https://phylipweb.github.io/phylip/download/${PHYLIP_TARBALL}"

echo "=== PHYLIP ${PHYLIP_VERSION} Setup ==="

# Check if already compiled
if [ -x "${PHYLIP_DIR}/exe/dnadist" ]; then
    echo "PHYLIP ${PHYLIP_VERSION} is already compiled."
    echo "Executables at: $(pwd)/${PHYLIP_DIR}/exe/"
    exit 0
fi

# Download
if [ ! -f "${PHYLIP_TARBALL}" ]; then
    echo "Downloading PHYLIP ${PHYLIP_VERSION}..."
    curl -L -o "${PHYLIP_TARBALL}" "${PHYLIP_URL}"
fi

# Extract
if [ ! -d "${PHYLIP_DIR}" ]; then
    echo "Extracting..."
    tar xzf "${PHYLIP_TARBALL}"
fi

# Compile
echo "Compiling PHYLIP ${PHYLIP_VERSION}..."
cd "${PHYLIP_DIR}/src"

# Use the Unix Makefile (works on macOS and Linux)
if [ -f "Makefile.unx" ]; then
    make -f Makefile.unx install
elif [ -f "Makefile.osx" ]; then
    make -f Makefile.osx install
else
    echo "ERROR: No suitable Makefile found"
    exit 1
fi

cd ../..
echo ""
echo "=== PHYLIP ${PHYLIP_VERSION} compiled successfully ==="
echo "Executables at: $(pwd)/${PHYLIP_DIR}/exe/"
echo ""
echo "Key programs for validation:"
ls -1 "${PHYLIP_DIR}/exe/" | grep -E '^(dnadist|dnaml|dnamlk|dnapars|neighbor|seqboot|consense|protdist|protpars)$' | while read prog; do
    echo "  ${prog}"
done
