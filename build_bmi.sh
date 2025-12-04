#!/bin/bash
#
# Build script for openLISEM BMI wrapper
# Usage: ./build_bmi.sh [clean|install|test]
#

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${SCRIPT_DIR}/build"

echo -e "${GREEN}=====================================${NC}"
echo -e "${GREEN}openLISEM BMI Wrapper Build Script${NC}"
echo -e "${GREEN}=====================================${NC}"

# Parse arguments
CLEAN=false
INSTALL=false
RUN_TESTS=false

for arg in "$@"; do
    case $arg in
        clean)
            CLEAN=true
            ;;
        install)
            INSTALL=true
            ;;
        test)
            RUN_TESTS=true
            ;;
        *)
            echo -e "${RED}Unknown argument: $arg${NC}"
            echo "Usage: $0 [clean|install|test]"
            exit 1
            ;;
    esac
done

# Clean if requested
if [ "$CLEAN" = true ]; then
    echo -e "${YELLOW}Cleaning build directory...${NC}"
    rm -rf "${BUILD_DIR}"
    echo -e "${GREEN}Clean complete!${NC}"
fi

# Check dependencies
echo -e "\n${YELLOW}Checking dependencies...${NC}"

command -v cmake >/dev/null 2>&1 || {
    echo -e "${RED}ERROR: cmake is required but not installed.${NC}"
    echo "Install with: sudo apt-get install cmake"
    exit 1
}

command -v python3 >/dev/null 2>&1 || {
    echo -e "${RED}ERROR: python3 is required but not installed.${NC}"
    exit 1
}

python3 -c "import pybind11" 2>/dev/null || {
    echo -e "${YELLOW}WARNING: pybind11 not found in Python.${NC}"
    echo "Install with: pip3 install pybind11"
    echo -e "${YELLOW}Continuing anyway...${NC}"
}

echo -e "${GREEN}Dependencies OK!${NC}"

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

# Configure
echo -e "\n${YELLOW}Configuring CMake...${NC}"
cmake .. -DBUILD_BMI=ON -DCMAKE_BUILD_TYPE=Release

# Build
echo -e "\n${YELLOW}Building...${NC}"
make -j$(nproc)

if [ $? -eq 0 ]; then
    echo -e "\n${GREEN}Build successful!${NC}"
    echo -e "BMI library location: ${BUILD_DIR}/bmi/libbmi_openlisem_core.so"
    echo -e "Python module location: ${BUILD_DIR}/bmi/bmi_openlisem*.so"
else
    echo -e "\n${RED}Build failed!${NC}"
    exit 1
fi

# Install if requested
if [ "$INSTALL" = true ]; then
    echo -e "\n${YELLOW}Installing Python module...${NC}"
    cd "${SCRIPT_DIR}"
    pip3 install -e .

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}Installation successful!${NC}"
    else
        echo -e "${RED}Installation failed!${NC}"
        exit 1
    fi
fi

# Run tests if requested
if [ "$RUN_TESTS" = true ]; then
    echo -e "\n${YELLOW}Running tests...${NC}"

    # Test Python import
    python3 -c "import bmi_openlisem; print('Python module import: OK')" || {
        echo -e "${RED}Python module import failed!${NC}"
        exit 1
    }

    # Test basic functionality
    python3 -c "
import bmi_openlisem
model = bmi_openlisem.BmiOpenLISEM()
print(f'Component name: {model.get_component_name()}')
print('Basic functionality: OK')
" || {
        echo -e "${RED}Basic functionality test failed!${NC}"
        exit 1
    }

    echo -e "${GREEN}All tests passed!${NC}"
fi

# Summary
echo -e "\n${GREEN}=====================================${NC}"
echo -e "${GREEN}Build Summary${NC}"
echo -e "${GREEN}=====================================${NC}"
echo -e "Build directory: ${BUILD_DIR}"
echo -e "BMI library: ${BUILD_DIR}/bmi/libbmi_openlisem_core.so"
echo -e "Python module: ${BUILD_DIR}/bmi/bmi_openlisem*.so"

echo -e "\n${YELLOW}To use the Python module:${NC}"
echo -e "  1. Install: pip3 install -e ."
echo -e "  2. Or add to PYTHONPATH: export PYTHONPATH=${BUILD_DIR}/bmi:\$PYTHONPATH"

echo -e "\n${YELLOW}To test:${NC}"
echo -e "  python3 examples/test_bmi_simple.py <runfile.run>"

echo -e "\n${GREEN}Done!${NC}"
