#!/usr/bin/env python3
"""
Advanced BMI test for model coupling experiments

This script demonstrates:
- Setting input variables (e.g., rainfall)
- Coupling with external models
- Advanced state manipulation
- Variable exchange between models
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

try:
    import bmi_openlisem
except ImportError:
    print("ERROR: bmi_openlisem module not found!")
    print("Please install it first using: pip install -e .")
    sys.exit(1)


def run_coupling_experiment(config_file):
    """
    Demonstrate model coupling capabilities

    This example shows how to:
    1. Initialize the model
    2. Modify rainfall inputs dynamically
    3. Exchange state variables
    4. Run coupled time stepping
    """

    print("=" * 60)
    print("openLISEM BMI - Model Coupling Example")
    print("=" * 60)

    # Initialize model
    model = bmi_openlisem.BmiOpenLISEM()
    model.initialize(config_file)

    grid_id = 0
    grid_shape = model.get_grid_shape(grid_id)
    nrows, ncols = grid_shape

    print(f"\nGrid: {nrows} x {ncols} cells")
    print(f"Time step: {model.get_time_step()} seconds")

    # Simulate coupling with another model
    # Example: Dynamic rainfall modification
    print("\n" + "=" * 60)
    print("Scenario: Coupling with rainfall model")
    print("=" * 60)

    times = []
    rainfall_values = []
    discharge_values = []
    infiltration_values = []

    # Run simulation with dynamic rainfall
    current_time = model.get_start_time()
    end_time = model.get_end_time()
    dt = model.get_time_step()

    step = 0
    while current_time < end_time:
        step += 1

        # Example: Modify rainfall based on external model
        # (In real coupling, this would come from another model)
        time_minutes = current_time / 60.0

        # Simulate a rainfall event that varies in time
        if time_minutes < 10:
            # Increasing rainfall
            rainfall_intensity = 5.0 + time_minutes * 2.0  # mm/h
        elif time_minutes < 30:
            # Peak rainfall
            rainfall_intensity = 25.0
        elif time_minutes < 45:
            # Decreasing rainfall
            rainfall_intensity = 25.0 - (time_minutes - 30) * 1.5
        else:
            # No rainfall
            rainfall_intensity = 0.0

        rainfall_values.append(rainfall_intensity)

        # Set rainfall (if the variable exists and is settable)
        try:
            rainfall_grid = np.ones((nrows, ncols)) * rainfall_intensity
            model.set_value("rainfall_intensity", rainfall_grid)
        except:
            pass  # Variable might not be settable in this implementation

        # Update model
        model.update()
        current_time = model.get_current_time()
        times.append(current_time)

        # Get state variables for analysis
        try:
            discharge = model.get_value("surface_runoff")
            infiltration = model.get_value("infiltration_rate")

            max_discharge = np.nanmax(discharge)
            mean_infiltration = np.nanmean(infiltration[infiltration > 0]) if np.any(infiltration > 0) else 0

            discharge_values.append(max_discharge)
            infiltration_values.append(mean_infiltration)

            if step % 10 == 0:
                print(f"Step {step:3d} | Time: {current_time/60:6.1f} min | "
                      f"Rainfall: {rainfall_intensity:5.1f} mm/h | "
                      f"Discharge: {max_discharge:6.3f} m³/s | "
                      f"Infiltration: {mean_infiltration:5.2f} mm/h")
        except Exception as e:
            if step % 10 == 0:
                print(f"Step {step:3d} | Time: {current_time/60:6.1f} min | "
                      f"Rainfall: {rainfall_intensity:5.1f} mm/h | "
                      f"Error: {e}")

        # Stop after a reasonable number of steps for testing
        if step >= 100:
            break

    # Plot coupling results
    times_min = np.array(times) / 60.0

    fig, axes = plt.subplots(3, 1, figsize=(12, 10))

    # Rainfall forcing
    axes[0].plot(times_min, rainfall_values, 'b-', linewidth=2)
    axes[0].set_ylabel('Rainfall (mm/h)')
    axes[0].set_title('Model Coupling: Dynamic Rainfall Input')
    axes[0].grid(True, alpha=0.3)

    # Discharge response
    if len(discharge_values) > 0:
        axes[1].plot(times_min, discharge_values, 'r-', linewidth=2)
        axes[1].set_ylabel('Discharge (m³/s)')
        axes[1].set_title('Model Response: Surface Runoff')
        axes[1].grid(True, alpha=0.3)

    # Infiltration response
    if len(infiltration_values) > 0:
        axes[2].plot(times_min, infiltration_values, 'g-', linewidth=2)
        axes[2].set_ylabel('Infiltration (mm/h)')
        axes[2].set_xlabel('Time (minutes)')
        axes[2].set_title('Model Response: Infiltration Rate')
        axes[2].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('bmi_coupling_results.png', dpi=150)
    print("\nCoupling results saved to: bmi_coupling_results.png")

    # Demonstrate direct variable access (zero-copy)
    print("\n" + "=" * 60)
    print("Demonstrating direct variable access (zero-copy):")
    print("=" * 60)

    try:
        # Get pointer to variable (no copy)
        water_ptr = model.get_value_ptr("surface_water_depth")
        print(f"Water depth array shape: {water_ptr.shape}")
        print(f"Water depth array dtype: {water_ptr.dtype}")
        print(f"Min water depth: {np.nanmin(water_ptr):.6f} m")
        print(f"Max water depth: {np.nanmax(water_ptr):.6f} m")
        print(f"Mean water depth: {np.nanmean(water_ptr[water_ptr > 0]):.6f} m")
    except Exception as e:
        print(f"Could not access water depth pointer: {e}")

    # Finalize
    model.finalize()
    print("\nModel coupling experiment completed!")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python test_bmi_coupling.py <runfile.run>")
        print("\nThis script demonstrates advanced BMI features for model coupling:")
        print("  - Dynamic input modification")
        print("  - State variable exchange")
        print("  - Zero-copy variable access")
        print("\nExample:")
        print("  python test_bmi_coupling.py /path/to/your/runfile.run")
        sys.exit(1)

    config_file = sys.argv[1]
    run_coupling_experiment(config_file)
