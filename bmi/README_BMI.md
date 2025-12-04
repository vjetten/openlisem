# openLISEM BMI Wrapper

This directory contains the Basic Model Interface (BMI) wrapper for openLISEM, enabling model coupling and integration with other modeling frameworks.

## Overview

The BMI (Basic Model Interface) is a standardized set of functions for model control and data exchange, developed by the CSDMS (Community Surface Dynamics Modeling System). This wrapper allows openLISEM to be used in:

- Model coupling frameworks (e.g., GLUE, eWaterCycle, SUMMA)
- Multi-model workflows
- Python-based modeling environments
- Automated calibration systems
- Cloud-based model services

## Features

- ✅ Full BMI 2.0 specification compliance
- ✅ All standard BMI functions (initialize, update, finalize, etc.)
- ✅ Variable getters/setters for key hydrological variables
- ✅ Grid information access
- ✅ Python bindings via pybind11
- ✅ C interface for maximum compatibility
- ✅ Zero-copy data access for performance
- ✅ Support for time stepping and state control

## Available Variables

### Input Variables
- `rainfall_intensity` - Rainfall intensity (mm/h)

### Output Variables
- `surface_water_depth` - Surface water depth (m)
- `surface_runoff` - Surface runoff discharge (m³/s)
- `surface_velocity` - Surface water velocity (m/s)
- `infiltration_rate` - Infiltration rate (mm/h)
- `cumulative_infiltration` - Cumulative infiltration (mm)
- `soil_moisture_content` - Soil moisture content (m³/m³)
- `detachment_rate` - Erosion detachment rate (kg/m²/s)
- `sediment_concentration` - Sediment concentration (kg/m³)
- `cumulative_erosion` - Cumulative erosion (kg/m²)
- `channel_discharge` - Channel discharge (m³/s)
- `channel_water_depth` - Channel water depth (m)

## Installation

### Prerequisites

- CMake >= 3.16
- Qt6 (Core, Gui, Network)
- GDAL
- OpenMP
- Python >= 3.7 (for Python bindings)
- pybind11 >= 2.6.0 (for Python bindings)
- NumPy >= 1.18.0 (for Python bindings)

### Building and Installing

#### Option 1: Python Installation (Recommended)

```bash
# From the repository root
pip install -e .
```

This will:
1. Compile the C++ BMI wrapper
2. Build the Python bindings
3. Install the `bmi_openlisem` Python module

#### Option 2: Manual CMake Build

```bash
# Create build directory
mkdir build && cd build

# Configure with BMI enabled
cmake .. -DBUILD_BMI=ON

# Build
make -j4

# Install (optional)
sudo make install
```

#### Option 3: Build Python module manually

```bash
# Install in development mode
python setup.py develop

# Or install normally
python setup.py install
```

## Usage Examples

### Python Usage

#### Basic Example

```python
import bmi_openlisem
import numpy as np

# Create model instance
model = bmi_openlisem.BmiOpenLISEM()

# Initialize with runfile
model.initialize("/path/to/runfile.run")

# Get model info
print(f"Model: {model.get_component_name()}")
print(f"Time step: {model.get_time_step()} {model.get_time_units()}")

# Run simulation
while model.get_current_time() < model.get_end_time():
    model.update()

    # Get water depth
    water_depth = model.get_value("surface_water_depth")
    print(f"Time: {model.get_current_time():.1f}, Max depth: {np.max(water_depth):.4f} m")

# Finalize
model.finalize()
```

#### Advanced Coupling Example

```python
import bmi_openlisem
import numpy as np

# Initialize two models for coupling
model1 = bmi_openlisem.BmiOpenLISEM()
model1.initialize("catchment1.run")

model2 = bmi_openlisem.BmiOpenLISEM()
model2.initialize("catchment2.run")

# Synchronize time steps
dt = min(model1.get_time_step(), model2.get_time_step())

while model1.get_current_time() < model1.get_end_time():
    # Update both models
    model1.update()
    model2.update()

    # Exchange variables between models
    # Example: outlet of model1 becomes inlet of model2
    discharge1 = model1.get_value("channel_discharge")

    # Process and set as input to model2
    # (implementation depends on coupling strategy)

    print(f"Time: {model1.get_current_time():.1f} s")

model1.finalize()
model2.finalize()
```

#### Zero-Copy Access

```python
# Get direct pointer to model data (no copy)
water_depth_ptr = model.get_value_ptr("surface_water_depth")

# Modify directly (changes affect model state)
water_depth_ptr[50, 50] = 0.5  # Set water depth at cell (50, 50)

# Changes are immediately reflected in next update
model.update()
```

### C++ Usage

```cpp
#include "bmi_openlisem.h"

int main() {
    bmi::BmiOpenLISEM model;

    // Initialize
    model.Initialize("runfile.run");

    // Get grid info
    int grid_shape[2];
    model.GetGridShape(0, grid_shape);

    // Allocate data
    std::vector<double> water_depth(grid_shape[0] * grid_shape[1]);

    // Run
    while (model.GetCurrentTime() < model.GetEndTime()) {
        model.Update();
        model.GetValue("surface_water_depth", water_depth.data());
    }

    // Finalize
    model.Finalize();

    return 0;
}
```

### C Interface Usage

```c
#include "bmi.h"

int main() {
    void* model = bmi_new();

    bmi_initialize(model, "runfile.run");

    double current_time;
    bmi_get_current_time(model, &current_time);

    bmi_update(model);

    bmi_finalize(model);
    bmi_delete(model);

    return 0;
}
```

## Testing

Example test scripts are provided in the `examples/` directory:

```bash
# Run basic BMI test
python examples/test_bmi_simple.py /path/to/runfile.run

# Run coupling experiment
python examples/test_bmi_coupling.py /path/to/runfile.run

# Run C++ test
cd build
./test_bmi /path/to/runfile.run
```

## BMI Function Reference

### Control Functions
- `initialize(config_file)` - Initialize model
- `update()` - Advance one time step
- `update_until(time)` - Advance to specific time
- `finalize()` - Clean up and finalize

### Information Functions
- `get_component_name()` - Get model name
- `get_input_var_names()` - List input variables
- `get_output_var_names()` - List output variables
- `get_var_type(name)` - Get variable data type
- `get_var_units(name)` - Get variable units
- `get_var_grid(name)` - Get variable grid ID

### Time Functions
- `get_start_time()` - Get simulation start time
- `get_end_time()` - Get simulation end time
- `get_current_time()` - Get current model time
- `get_time_step()` - Get time step size
- `get_time_units()` - Get time units

### Variable Access
- `get_value(name)` - Get variable values (copy)
- `get_value_ptr(name)` - Get variable pointer (no copy)
- `set_value(name, values)` - Set variable values
- `get_value_at_indices(name, indices)` - Get values at specific indices
- `set_value_at_indices(name, indices, values)` - Set values at specific indices

### Grid Functions
- `get_grid_rank(grid_id)` - Get grid dimensions
- `get_grid_size(grid_id)` - Get total grid cells
- `get_grid_type(grid_id)` - Get grid type
- `get_grid_shape(grid_id)` - Get grid shape
- `get_grid_spacing(grid_id)` - Get cell size
- `get_grid_origin(grid_id)` - Get grid origin

## Architecture

```
bmi/
├── bmi.h                  # C BMI header
├── bmi_openlisem.h        # C++ BMI class header
├── bmi_openlisem.cpp      # C++ BMI implementation
├── bmi_python.cpp         # Python bindings (pybind11)
├── CMakeLists.txt         # Build configuration
└── README_BMI.md          # This file
```

## Known Limitations

1. **Headless Mode Only**: The BMI wrapper runs openLISEM in headless mode (no GUI)
2. **Single Timestep Control**: The model is designed for full simulation runs, so single-timestep control is implemented by breaking up the time loop
3. **Variable Availability**: Not all internal model variables are exposed through BMI (can be extended as needed)
4. **Configuration**: Model configuration is still done via runfiles (.run), not through BMI setters

## Extending the Wrapper

To add new variables to the BMI interface:

1. Register the variable in `InitializeVariableRegistry()`:
```cpp
RegisterVariable("new_variable", "double", "units", "node", 0, false, true);
```

2. Add pointer mapping in `GetVariablePointer()`:
```cpp
else if (name == "new_variable" && model_->NewVar) {
    return model_->NewVar->data[0];
}
```

## References

- BMI Specification: https://bmi.readthedocs.io/
- CSDMS: https://csdms.colorado.edu/
- openLISEM: https://github.com/vjetten/openlisem
- pybind11: https://pybind11.readthedocs.io/

## License

This BMI wrapper follows the same GPLv3 license as openLISEM.

## Support

For issues and questions:
- openLISEM issues: https://github.com/vjetten/openlisem/issues
- Contact: v.g.jetten@utwente.nl
