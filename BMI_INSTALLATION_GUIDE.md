# openLISEM BMI Wrapper - Installation Guide

This guide will help you build and install the BMI (Basic Model Interface) wrapper for openLISEM.

## Prerequisites

### System Requirements

- Linux (tested on Ubuntu 20.04+) or macOS
- GCC/Clang with C++17 support
- CMake >= 3.16
- Python >= 3.7 (for Python bindings)

### Required Dependencies

#### Ubuntu/Debian
```bash
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    cmake \
    qt6-base-dev \
    libgdal-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libomp-dev \
    python3-dev \
    python3-pip \
    python3-numpy \
    git
```

#### Install pybind11
```bash
pip3 install pybind11
```

Or install from source:
```bash
git clone https://github.com/pybind/pybind11.git
cd pybind11
mkdir build && cd build
cmake ..
make install
```

#### QWT (optional, for GUI version)
The BMI wrapper runs in headless mode and doesn't require QWT, but if you want to build the full GUI version:

```bash
# Download and build QWT with multi-axes support
cd /tmp
git clone https://github.com/osakared/qwt.git qwt-multiaxes
cd qwt-multiaxes
mkdir build && cd build
cmake ..
make -j4
sudo make install
```

## Installation Methods

### Method 1: Quick Python Installation (Recommended)

This is the easiest method and automatically handles the build process:

```bash
# Clone the repository (if you haven't already)
git clone https://github.com/vjetten/openlisem.git
cd openlisem

# Install in editable mode
pip3 install -e .
```

This will:
1. Compile the C++ BMI wrapper
2. Build the Python bindings
3. Install the `bmi_openlisem` Python module in your Python environment

### Method 2: Manual CMake Build

If you want more control over the build process:

```bash
# Navigate to the repository root
cd /path/to/openlisem_bmi

# Create build directory
mkdir build
cd build

# Configure with BMI enabled
cmake .. -DBUILD_BMI=ON

# Build
make -j$(nproc)

# The BMI library will be in: build/bmi/libbmi_openlisem_core.so
# The Python module will be in: build/bmi/bmi_openlisem.so
```

To install system-wide:
```bash
sudo make install
```

### Method 3: Build Python Module Manually

```bash
# From repository root
python3 setup.py build

# Install
python3 setup.py install

# Or install in development mode (changes reflected immediately)
python3 setup.py develop
```

## Verifying Installation

### Test Python Module

```python
import bmi_openlisem
print(bmi_openlisem.__file__)

# Create model instance
model = bmi_openlisem.BmiOpenLISEM()
print(f"Model component: {model.get_component_name()}")
```

### Run Example Scripts

```bash
# Make sure you have a valid runfile
# Download example data if needed from openLISEM repository

# Run basic test
python3 examples/test_bmi_simple.py /path/to/your/runfile.run

# Run coupling test
python3 examples/test_bmi_coupling.py /path/to/your/runfile.run
```

### Test C++ Library

```bash
cd build
g++ -o test_bmi ../examples/test_bmi.cpp \
    -I../bmi -I../include \
    -L./bmi -lbmi_openlisem_core \
    -lQt6Core -lgdal -fopenmp \
    -Wl,-rpath,./bmi

./test_bmi /path/to/your/runfile.run
```

## Troubleshooting

### Issue: CMake can't find Qt6

```bash
# Set Qt6 path
export Qt6_DIR=/path/to/qt6/lib/cmake/Qt6

# Or specify during configuration
cmake .. -DBUILD_BMI=ON -DQt6_DIR=/path/to/qt6/lib/cmake/Qt6
```

### Issue: Can't find pybind11

```bash
# Install via pip
pip3 install pybind11

# Or specify path
cmake .. -DBUILD_BMI=ON -Dpybind11_DIR=/path/to/pybind11/share/cmake/pybind11
```

### Issue: GDAL not found

```bash
# Ubuntu/Debian
sudo apt-get install libgdal-dev

# macOS
brew install gdal

# Or specify path
cmake .. -DBUILD_BMI=ON -DGDAL_DIR=/path/to/gdal
```

### Issue: Import error in Python

```bash
# Check if the module is in Python path
python3 -c "import sys; print('\n'.join(sys.path))"

# Add to PYTHONPATH if needed
export PYTHONPATH=/path/to/build/bmi:$PYTHONPATH

# Or reinstall
pip3 install -e . --force-reinstall
```

### Issue: Runtime library errors

```bash
# Make sure all shared libraries are in library path
export LD_LIBRARY_PATH=/path/to/build/bmi:$LD_LIBRARY_PATH

# Check dependencies
ldd /path/to/bmi_openlisem.so
```

### Issue: Model initialization fails

Common causes:
1. **Invalid runfile path**: Make sure the path to the .run file is correct
2. **Missing input files**: Check that all maps and input files specified in the runfile exist
3. **Permission issues**: Ensure you have read/write permissions for input/output directories
4. **GDAL configuration**: Some systems need `GDAL_DATA` environment variable set

```bash
export GDAL_DATA=/usr/share/gdal
```

## Building for Development

If you're developing the BMI wrapper:

```bash
# Build in debug mode
cmake .. -DBUILD_BMI=ON -DCMAKE_BUILD_TYPE=Debug

# Enable verbose output
make VERBOSE=1

# Run with debugging
gdb --args python3 examples/test_bmi_simple.py runfile.run
```

## Docker Installation (Alternative)

Create a Dockerfile:

```dockerfile
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    build-essential cmake git \
    qt6-base-dev libgdal-dev \
    libssl-dev libcurl4-openssl-dev libomp-dev \
    python3-dev python3-pip python3-numpy

RUN pip3 install pybind11

WORKDIR /app
COPY . .

RUN pip3 install -e .

CMD ["python3"]
```

Build and run:
```bash
docker build -t openlisem-bmi .
docker run -it -v /path/to/data:/data openlisem-bmi
```

## Performance Tips

1. **OpenMP Threads**: Control parallelization
   ```bash
   export OMP_NUM_THREADS=4
   ```

2. **Build Optimization**: Use Release mode
   ```bash
   cmake .. -DBUILD_BMI=ON -DCMAKE_BUILD_TYPE=Release
   ```

3. **Link-Time Optimization**:
   ```bash
   cmake .. -DBUILD_BMI=ON -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
   ```

## Next Steps

After successful installation:

1. Review the [BMI README](bmi/README_BMI.md) for detailed API documentation
2. Try the example scripts in the `examples/` directory
3. Check the available variables with:
   ```python
   model = bmi_openlisem.BmiOpenLISEM()
   model.initialize("runfile.run")
   print(model.get_output_var_names())
   ```
4. Start building your model coupling experiments!

## Getting Help

- openLISEM GitHub: https://github.com/vjetten/openlisem
- BMI Documentation: https://bmi.readthedocs.io/
- Email: v.g.jetten@utwente.nl

## Known Limitations

1. The BMI wrapper requires openLISEM to run in headless mode (no GUI)
2. Some advanced GUI features are not accessible through BMI
3. Model configuration must be provided via runfile, not through BMI functions
4. Thread safety: Create separate model instances for parallel runs

## Contributing

To contribute improvements to the BMI wrapper:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

Key areas for contribution:
- Adding more exposed variables
- Improving error handling
- Adding unit tests
- Performance optimizations
- Documentation improvements
