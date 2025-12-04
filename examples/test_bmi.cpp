/*************************************************************************
**  C++ example for using openLISEM BMI wrapper
**
**  Compile with:
**    g++ -o test_bmi test_bmi.cpp -I../bmi -L../build/bmi \
**        -lbmi_openlisem_core -lQt6Core -lgdal -fopenmp
**************************************************************************/

#include <iostream>
#include <vector>
#include <cmath>
#include "bmi_openlisem.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <runfile.run>" << std::endl;
        return 1;
    }

    std::string config_file = argv[1];

    std::cout << "========================================" << std::endl;
    std::cout << "openLISEM BMI C++ Test" << std::endl;
    std::cout << "========================================" << std::endl;

    try {
        // Create model instance
        bmi::BmiOpenLISEM model;

        // Initialize
        std::cout << "\nInitializing model with: " << config_file << std::endl;
        if (model.Initialize(config_file) != BMI_SUCCESS) {
            std::cerr << "Failed to initialize model!" << std::endl;
            return 1;
        }

        // Get model information
        std::cout << "\nModel: " << model.GetComponentName() << std::endl;
        std::cout << "Start time: " << model.GetStartTime() << " "
                  << model.GetTimeUnits() << std::endl;
        std::cout << "End time: " << model.GetEndTime() << " "
                  << model.GetTimeUnits() << std::endl;
        std::cout << "Time step: " << model.GetTimeStep() << " "
                  << model.GetTimeUnits() << std::endl;

        // Get grid information
        int grid_id = 0;
        int grid_shape[2];
        double grid_spacing[2];
        double grid_origin[2];

        model.GetGridShape(grid_id, grid_shape);
        model.GetGridSpacing(grid_id, grid_spacing);
        model.GetGridOrigin(grid_id, grid_origin);

        std::cout << "\nGrid information:" << std::endl;
        std::cout << "  Shape: " << grid_shape[0] << " x " << grid_shape[1] << std::endl;
        std::cout << "  Spacing: " << grid_spacing[0] << " m" << std::endl;
        std::cout << "  Origin: (" << grid_origin[0] << ", " << grid_origin[1] << ")" << std::endl;

        // Get available variables
        auto input_vars = model.GetInputVarNames();
        auto output_vars = model.GetOutputVarNames();

        std::cout << "\nInput variables (" << input_vars.size() << "):" << std::endl;
        for (const auto& var : input_vars) {
            std::cout << "  - " << var << " [" << model.GetVarUnits(var) << "]" << std::endl;
        }

        std::cout << "\nOutput variables (" << output_vars.size() << "):" << std::endl;
        for (const auto& var : output_vars) {
            std::cout << "  - " << var << " [" << model.GetVarUnits(var) << "]" << std::endl;
        }

        // Run the model
        std::cout << "\n========================================" << std::endl;
        std::cout << "Running model..." << std::endl;
        std::cout << "========================================" << std::endl;

        int n_steps = 10;
        int grid_size = grid_shape[0] * grid_shape[1];
        std::vector<double> water_depth(grid_size);

        for (int step = 0; step < n_steps; step++) {
            // Update model
            if (model.Update() != BMI_SUCCESS) {
                std::cerr << "Error updating model at step " << step << std::endl;
                break;
            }

            double current_time = model.GetCurrentTime();

            // Get water depth
            try {
                model.GetValue("surface_water_depth", water_depth.data());

                // Calculate statistics
                double max_depth = 0.0;
                double sum_depth = 0.0;
                int count = 0;

                for (double depth : water_depth) {
                    if (!std::isnan(depth) && depth > 0) {
                        max_depth = std::max(max_depth, depth);
                        sum_depth += depth;
                        count++;
                    }
                }

                double mean_depth = count > 0 ? sum_depth / count : 0.0;

                std::cout << "Step " << (step + 1) << " | Time: " << current_time
                          << " s | Max depth: " << max_depth
                          << " m | Mean depth: " << mean_depth << " m" << std::endl;
            }
            catch (const std::exception& e) {
                std::cout << "Step " << (step + 1) << " | Time: " << current_time
                          << " s | Error: " << e.what() << std::endl;
            }

            if (current_time >= model.GetEndTime()) {
                std::cout << "\nReached end time!" << std::endl;
                break;
            }
        }

        // Finalize
        std::cout << "\n========================================" << std::endl;
        std::cout << "Finalizing model..." << std::endl;
        model.Finalize();
        std::cout << "Done!" << std::endl;
        std::cout << "========================================" << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
