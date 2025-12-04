#!/usr/bin/env python3
"""
Simple test script for openLISEM BMI wrapper

This script demonstrates basic BMI functionality:
- Initialize the model
- Step through time
- Access model variables
- Finalize the model
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


def run_bmi_example(config_file):
    """
    Run a simple BMI example with openLISEM

    Parameters
    ----------
    config_file : str
        Path to the openLISEM runfile (.run)
    """

    print("=" * 60)
    print("openLISEM BMI Test")
    print("=" * 60)

    # Create model instance
    model = bmi_openlisem.BmiOpenLISEM()

    # Initialize
    print(f"\nInitializing model with: {config_file}")
    model.initialize(config_file)

    # Get model information
    print(f"\nModel component: {model.get_component_name()}")
    print(f"Start time: {model.get_start_time()} {model.get_time_units()}")
    print(f"End time: {model.get_end_time()} {model.get_time_units()}")
    print(f"Time step: {model.get_time_step()} {model.get_time_units()}")

    # Get grid information
    grid_id = 0
    grid_shape = model.get_grid_shape(grid_id)
    grid_spacing = model.get_grid_spacing(grid_id)
    grid_origin = model.get_grid_origin(grid_id)

    print(f"\nGrid information:")
    print(f"  Shape: {grid_shape} (rows x cols)")
    print(f"  Spacing: {grid_spacing} m")
    print(f"  Origin: {grid_origin}")
    print(f"  Total cells: {model.get_grid_size(grid_id)}")

    # Get available variables
    print(f"\nInput variables ({model.get_input_item_count()}):")
    for var in model.get_input_var_names():
        print(f"  - {var} [{model.get_var_units(var)}]")

    print(f"\nOutput variables ({model.get_output_item_count()}):")
    for var in model.get_output_var_names():
        print(f"  - {var} [{model.get_var_units(var)}]")

    # Run the model
    print("\n" + "=" * 60)
    print("Running model...")
    print("=" * 60)

    n_steps = 10  # Run for 10 time steps
    times = []
    outlet_discharge = []

    for step in range(n_steps):
        # Update model by one time step
        model.update()

        current_time = model.get_current_time()
        times.append(current_time)

        # Get surface water depth
        try:
            water_depth = model.get_value("surface_water_depth")
            max_depth = np.nanmax(water_depth)
            mean_depth = np.nanmean(water_depth[water_depth > 0]) if np.any(water_depth > 0) else 0

            # Get discharge
            discharge = model.get_value("surface_runoff")
            max_discharge = np.nanmax(discharge)
            outlet_discharge.append(max_discharge)

            print(f"Step {step+1:3d} | Time: {current_time:8.1f} s | "
                  f"Max water depth: {max_depth:.4f} m | "
                  f"Mean water depth: {mean_depth:.4f} m | "
                  f"Max discharge: {max_discharge:.4f} m³/s")
        except Exception as e:
            print(f"Step {step+1:3d} | Time: {current_time:8.1f} s | Error accessing variables: {e}")

        # Stop if we reached the end
        if current_time >= model.get_end_time():
            print("\nReached end time!")
            break

    # Plot results
    if len(times) > 0 and len(outlet_discharge) > 0:
        plt.figure(figsize=(10, 6))
        plt.plot(np.array(times) / 60, outlet_discharge, 'b-', linewidth=2)
        plt.xlabel('Time (minutes)')
        plt.ylabel('Maximum Discharge (m³/s)')
        plt.title('openLISEM BMI - Discharge Hydrograph')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig('bmi_discharge_hydrograph.png', dpi=150)
        print("\nHydrograph saved to: bmi_discharge_hydrograph.png")

    # Get final state
    print("\n" + "=" * 60)
    print("Final model state:")
    print("=" * 60)

    try:
        water_depth = model.get_value("surface_water_depth")
        infiltration = model.get_value("cumulative_infiltration")

        print(f"Total water on surface: {np.nansum(water_depth):.2f} m³")
        print(f"Mean infiltration: {np.nanmean(infiltration[infiltration > 0]):.2f} mm")

        # Plot final water depth
        plt.figure(figsize=(12, 8))
        plt.subplot(1, 2, 1)
        plt.imshow(water_depth, cmap='Blues', interpolation='nearest')
        plt.colorbar(label='Water depth (m)')
        plt.title('Final Surface Water Depth')

        plt.subplot(1, 2, 2)
        plt.imshow(infiltration, cmap='YlGn', interpolation='nearest')
        plt.colorbar(label='Infiltration (mm)')
        plt.title('Cumulative Infiltration')

        plt.tight_layout()
        plt.savefig('bmi_final_state.png', dpi=150)
        print("Final state maps saved to: bmi_final_state.png")
    except Exception as e:
        print(f"Error creating final state plots: {e}")

    # Finalize
    print("\n" + "=" * 60)
    print("Finalizing model...")
    model.finalize()
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python test_bmi_simple.py <runfile.run>")
        print("\nExample:")
        print("  python test_bmi_simple.py /path/to/your/runfile.run")
        sys.exit(1)

    config_file = sys.argv[1]
    run_bmi_example(config_file)
