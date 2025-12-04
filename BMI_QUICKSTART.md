# openLISEM BMI Wrapper - Quick Start Guide

Get started with the openLISEM BMI wrapper in 5 minutes!

## 🚀 Quick Installation

```bash
# Install dependencies (Ubuntu/Debian)
sudo apt-get install -y build-essential cmake qt6-base-dev libgdal-dev \
    libssl-dev libcurl4-openssl-dev libomp-dev python3-dev python3-pip python3-numpy

# Install pybind11
pip3 install pybind11

# Install openLISEM BMI wrapper
pip3 install -e .
```

Or use the build script:
```bash
./build_bmi.sh install
```

## 🎯 Basic Usage

### Python Example (5 lines!)

```python
import bmi_openlisem

model = bmi_openlisem.BmiOpenLISEM()
model.initialize("path/to/runfile.run")
model.update()  # Run one timestep
water_depth = model.get_value("surface_water_depth")
model.finalize()
```

### Full Simulation

```python
import bmi_openlisem
import numpy as np

# Initialize
model = bmi_openlisem.BmiOpenLISEM()
model.initialize("runfile.run")

# Run simulation
while model.get_current_time() < model.get_end_time():
    model.update()

    # Access results
    wh = model.get_value("surface_water_depth")
    q = model.get_value("surface_runoff")

    print(f"Time: {model.get_current_time():.1f} s, "
          f"Max water depth: {np.max(wh):.3f} m")

# Finalize
model.finalize()
```

## 📊 Available Variables

### Outputs (Read)
- `surface_water_depth` (m)
- `surface_runoff` (m³/s)
- `infiltration_rate` (mm/h)
- `cumulative_infiltration` (mm)
- `soil_moisture_content` (m³/m³)
- `sediment_concentration` (kg/m³)
- `channel_discharge` (m³/s)
- And more...

### Inputs (Write)
- `rainfall_intensity` (mm/h)

## 🔧 Common Operations

### Get Model Info
```python
print(model.get_component_name())          # "openLISEM"
print(model.get_time_step())               # e.g., 60.0 (seconds)
print(model.get_output_var_names())        # List all variables
```

### Grid Information
```python
grid_id = 0
shape = model.get_grid_shape(grid_id)      # (rows, cols)
spacing = model.get_grid_spacing(grid_id)  # cell size
origin = model.get_grid_origin(grid_id)    # (x, y)
```

### Set Values
```python
# Modify rainfall
rainfall = np.ones(grid_shape) * 25.0  # 25 mm/h
model.set_value("rainfall_intensity", rainfall)
```

### Zero-Copy Access (Fast!)
```python
# Get direct pointer (no copy)
water_ptr = model.get_value_ptr("surface_water_depth")

# Modify directly
water_ptr[50, 50] = 1.0  # Set water at cell (50,50)
```

## 🧪 Test Installation

```bash
# Run example script
python3 examples/test_bmi_simple.py examples/example.run

# Or test import
python3 -c "import bmi_openlisem; print('OK!')"
```

## 🔗 Model Coupling Example

```python
import bmi_openlisem

# Two coupled catchments
model1 = bmi_openlisem.BmiOpenLISEM()
model2 = bmi_openlisem.BmiOpenLISEM()

model1.initialize("catchment1.run")
model2.initialize("catchment2.run")

# Coupled time stepping
while model1.get_current_time() < model1.get_end_time():
    # Update both
    model1.update()
    model2.update()

    # Exchange data
    outflow1 = model1.get_value("channel_discharge")
    # ... process and pass to model2 ...

    print(f"Time: {model1.get_current_time()}")

model1.finalize()
model2.finalize()
```

## 📁 Example Files

- `examples/test_bmi_simple.py` - Basic usage
- `examples/test_bmi_coupling.py` - Advanced coupling
- `examples/test_bmi.cpp` - C++ example

## 🐛 Troubleshooting

### Can't import module
```bash
pip3 install -e . --force-reinstall
```

### Library not found
```bash
export LD_LIBRARY_PATH=./build/bmi:$LD_LIBRARY_PATH
```

### GDAL errors
```bash
export GDAL_DATA=/usr/share/gdal
```

## 📚 Documentation

- Full API: `bmi/README_BMI.md`
- Installation Guide: `BMI_INSTALLATION_GUIDE.md`
- BMI Specification: https://bmi.readthedocs.io/

## 💡 Tips

1. **Use `get_value_ptr()` for large grids** - avoids copying data
2. **Check variable names** - use `get_output_var_names()`
3. **Time units are seconds** - convert to minutes/hours as needed
4. **Always call `finalize()`** - cleans up resources
5. **Run in headless mode** - no GUI required

## 🎓 Learn More

See the full examples directory for:
- Rainfall-runoff simulations
- Erosion modeling
- Multi-model coupling
- Time series analysis
- Spatial output visualization

## 🆘 Help

- GitHub: https://github.com/vjetten/openlisem
- Email: v.g.jetten@utwente.nl
- BMI Community: https://csdms.colorado.edu/

---

**Ready to go? Run this now:**

```bash
pip3 install -e .
python3 -c "import bmi_openlisem; print('BMI wrapper installed successfully!')"
```
