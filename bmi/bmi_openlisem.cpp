/*************************************************************************
**  BMI wrapper implementation for openLISEM
**  Copyright (C) 2024
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**************************************************************************/

#include "bmi_openlisem.h"
#include "model.h"
#include "global.h"
#include "fixture.h"
#include <QCoreApplication>
#include <stdexcept>
#include <cstring>
#include <iostream>

namespace bmi {

BmiOpenLISEM::BmiOpenLISEM()
    : model_(nullptr),
      qapp_(nullptr),
      current_time_(0.0),
      start_time_(0.0),
      end_time_(0.0),
      time_step_(0.0),
      initialized_(false),
      finalized_(false),
      nrows_(0),
      ncols_(0),
      cellsize_(0.0),
      xll_(0.0),
      yll_(0.0) {
}

BmiOpenLISEM::~BmiOpenLISEM() {
    if (initialized_ && !finalized_) {
        Finalize();
    }
}

int BmiOpenLISEM::Initialize(const std::string& config_file) {
    try {
        if (initialized_) {
            std::cerr << "Model already initialized!" << std::endl;
            return BMI_FAILURE;
        }

        config_file_ = config_file;

        // Initialize Qt application (required for openLISEM)
        int argc = 1;
        char* argv[] = {(char*)"openlisem_bmi"};
        qapp_ = std::make_unique<QCoreApplication>(argc, argv);

        // Initialize GDAL
        static Fixture gdal_fixture;

        // Create model instance
        model_ = std::make_unique<TWorld>();

        // Set headless mode
        model_->noInterface = true;

        // Read the runfile
        QString qConfigFile = QString::fromStdString(config_file);
        if (!model_->readRunfile(qConfigFile)) {
            throw std::runtime_error("Failed to read runfile: " + config_file);
        }

        // Initialize the model (calls ReadMapData and other initialization)
        if (!model_->runModelSetup()) {
            throw std::runtime_error("Failed to initialize model");
        }

        // Get time information from model
        start_time_ = model_->BeginTime;
        end_time_ = model_->EndTime;
        time_step_ = model_->_dt;
        current_time_ = start_time_;

        // Get grid information
        nrows_ = model_->_nrRows;
        ncols_ = model_->_nrCols;
        cellsize_ = model_->_dx;

        // Get origin from DEM if available
        if (model_->DEM) {
            xll_ = model_->DEM->west();
            yll_ = model_->DEM->south();
        }

        // Initialize variable registry
        InitializeVariableRegistry();

        initialized_ = true;
        std::cout << "openLISEM initialized successfully via BMI" << std::endl;
        std::cout << "  Start time: " << start_time_ << " s" << std::endl;
        std::cout << "  End time: " << end_time_ << " s" << std::endl;
        std::cout << "  Time step: " << time_step_ << " s" << std::endl;
        std::cout << "  Grid: " << nrows_ << " x " << ncols_ << " cells" << std::endl;

        return BMI_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "Error in Initialize: " << e.what() << std::endl;
        return BMI_FAILURE;
    }
}

int BmiOpenLISEM::Update() {
    try {
        if (!initialized_ || finalized_) {
            throw std::runtime_error("Model not initialized or already finalized");
        }

        if (current_time_ >= end_time_) {
            std::cout << "Model has reached end time" << std::endl;
            return BMI_SUCCESS;
        }

        // Run one timestep
        // Note: openLISEM's DoModel() runs the entire simulation
        // For BMI, we need to step through one timestep at a time
        // This requires modifying the time loop to be externally controllable

        // Store the original time values
        double orig_begin_time = model_->BeginTime;
        double orig_end_time = model_->EndTime;

        // Set the model to run for just one timestep
        model_->BeginTime = current_time_;
        model_->EndTime = current_time_ + time_step_;

        // Execute one model step
        model_->time = current_time_;

        // Run the core model processes for this timestep
        model_->GetInputTimeseries();
        model_->InfilDynamicCrusting();
        model_->HydrologyProcesses();
        model_->ToTiledrain();
        model_->OverlandFlow();
        model_->ChannelFlowandErosion();
        model_->TileFlow();
        model_->TotalsHydro();
        model_->TotalsFlow();
        model_->TotalsSediment();
        model_->MassBalance();

        // Update current time
        current_time_ += time_step_;
        model_->time = current_time_;

        // Restore original time bounds
        model_->BeginTime = orig_begin_time;
        model_->EndTime = orig_end_time;

        // Update variable pointers in case they changed
        UpdateVariablePointers();

        return BMI_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "Error in Update: " << e.what() << std::endl;
        return BMI_FAILURE;
    }
}

int BmiOpenLISEM::UpdateUntil(double time) {
    try {
        if (!initialized_ || finalized_) {
            throw std::runtime_error("Model not initialized or already finalized");
        }

        if (time < current_time_) {
            throw std::runtime_error("Target time is before current time");
        }

        if (time > end_time_) {
            throw std::runtime_error("Target time is beyond end time");
        }

        while (current_time_ < time) {
            int status = Update();
            if (status != BMI_SUCCESS) {
                return status;
            }
        }

        return BMI_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "Error in UpdateUntil: " << e.what() << std::endl;
        return BMI_FAILURE;
    }
}

int BmiOpenLISEM::Finalize() {
    try {
        if (!initialized_) {
            std::cerr << "Model not initialized" << std::endl;
            return BMI_FAILURE;
        }

        if (finalized_) {
            std::cerr << "Model already finalized" << std::endl;
            return BMI_FAILURE;
        }

        // Write final outputs if needed
        if (model_) {
            model_->reportToFile();
        }

        // Clean up
        model_.reset();
        qapp_.reset();

        finalized_ = true;
        std::cout << "openLISEM finalized successfully" << std::endl;

        return BMI_SUCCESS;
    }
    catch (const std::exception& e) {
        std::cerr << "Error in Finalize: " << e.what() << std::endl;
        return BMI_FAILURE;
    }
}

std::string BmiOpenLISEM::GetComponentName() const {
    return "openLISEM";
}

int BmiOpenLISEM::GetInputItemCount() const {
    return static_cast<int>(input_var_names_.size());
}

int BmiOpenLISEM::GetOutputItemCount() const {
    return static_cast<int>(output_var_names_.size());
}

std::vector<std::string> BmiOpenLISEM::GetInputVarNames() const {
    return input_var_names_;
}

std::vector<std::string> BmiOpenLISEM::GetOutputVarNames() const {
    return output_var_names_;
}

int BmiOpenLISEM::GetVarGrid(const std::string& name) const {
    auto it = var_info_.find(name);
    if (it != var_info_.end()) {
        return it->second.grid_id;
    }
    return -1;
}

std::string BmiOpenLISEM::GetVarType(const std::string& name) const {
    auto it = var_info_.find(name);
    if (it != var_info_.end()) {
        return it->second.type;
    }
    return "unknown";
}

std::string BmiOpenLISEM::GetVarUnits(const std::string& name) const {
    auto it = var_info_.find(name);
    if (it != var_info_.end()) {
        return it->second.units;
    }
    return "unknown";
}

int BmiOpenLISEM::GetVarItemsize(const std::string& name) const {
    std::string type = GetVarType(name);
    if (type == "double") return sizeof(double);
    if (type == "float") return sizeof(float);
    if (type == "int") return sizeof(int);
    return 0;
}

int BmiOpenLISEM::GetVarNbytes(const std::string& name) const {
    int grid = GetVarGrid(name);
    if (grid == 0) {
        return GetVarItemsize(name) * nrows_ * ncols_;
    }
    return 0;
}

std::string BmiOpenLISEM::GetVarLocation(const std::string& name) const {
    auto it = var_info_.find(name);
    if (it != var_info_.end()) {
        return it->second.location;
    }
    return "unknown";
}

double BmiOpenLISEM::GetCurrentTime() const {
    return current_time_;
}

double BmiOpenLISEM::GetStartTime() const {
    return start_time_;
}

double BmiOpenLISEM::GetEndTime() const {
    return end_time_;
}

double BmiOpenLISEM::GetTimeStep() const {
    return time_step_;
}

std::string BmiOpenLISEM::GetTimeUnits() const {
    return "seconds";
}

void BmiOpenLISEM::GetValue(const std::string& name, void* dest) {
    if (!initialized_ || finalized_) {
        throw std::runtime_error("Model not initialized or already finalized");
    }

    void* src = GetVariablePointer(name);
    if (!src) {
        throw std::runtime_error("Variable not found: " + name);
    }

    int nbytes = GetVarNbytes(name);
    std::memcpy(dest, src, nbytes);
}

void* BmiOpenLISEM::GetValuePtr(const std::string& name) {
    if (!initialized_ || finalized_) {
        throw std::runtime_error("Model not initialized or already finalized");
    }

    return GetVariablePointer(name);
}

void BmiOpenLISEM::GetValueAtIndices(const std::string& name, void* dest,
                                      int* inds, int count) {
    if (!initialized_ || finalized_) {
        throw std::runtime_error("Model not initialized or already finalized");
    }

    void* src = GetVariablePointer(name);
    if (!src) {
        throw std::runtime_error("Variable not found: " + name);
    }

    int itemsize = GetVarItemsize(name);
    char* src_ptr = static_cast<char*>(src);
    char* dest_ptr = static_cast<char*>(dest);

    for (int i = 0; i < count; i++) {
        std::memcpy(dest_ptr + i * itemsize,
                   src_ptr + inds[i] * itemsize,
                   itemsize);
    }
}

void BmiOpenLISEM::SetValue(const std::string& name, void* src) {
    if (!initialized_ || finalized_) {
        throw std::runtime_error("Model not initialized or already finalized");
    }

    void* dest = GetVariablePointer(name);
    if (!dest) {
        throw std::runtime_error("Variable not found: " + name);
    }

    int nbytes = GetVarNbytes(name);
    std::memcpy(dest, src, nbytes);
}

void BmiOpenLISEM::SetValueAtIndices(const std::string& name, int* inds,
                                      int count, void* src) {
    if (!initialized_ || finalized_) {
        throw std::runtime_error("Model not initialized or already finalized");
    }

    void* dest = GetVariablePointer(name);
    if (!dest) {
        throw std::runtime_error("Variable not found: " + name);
    }

    int itemsize = GetVarItemsize(name);
    char* dest_ptr = static_cast<char*>(dest);
    char* src_ptr = static_cast<char*>(src);

    for (int i = 0; i < count; i++) {
        std::memcpy(dest_ptr + inds[i] * itemsize,
                   src_ptr + i * itemsize,
                   itemsize);
    }
}

int BmiOpenLISEM::GetGridRank(int grid) const {
    if (grid == 0) return 2;  // 2D grid
    return -1;
}

int BmiOpenLISEM::GetGridSize(int grid) const {
    if (grid == 0) return nrows_ * ncols_;
    return -1;
}

std::string BmiOpenLISEM::GetGridType(int grid) const {
    if (grid == 0) return "uniform_rectilinear";
    return "unknown";
}

void BmiOpenLISEM::GetGridShape(int grid, int* shape) const {
    if (grid == 0) {
        shape[0] = nrows_;
        shape[1] = ncols_;
    }
}

void BmiOpenLISEM::GetGridSpacing(int grid, double* spacing) const {
    if (grid == 0) {
        spacing[0] = cellsize_;
        spacing[1] = cellsize_;
    }
}

void BmiOpenLISEM::GetGridOrigin(int grid, double* origin) const {
    if (grid == 0) {
        origin[0] = yll_;
        origin[1] = xll_;
    }
}

void BmiOpenLISEM::GetGridX(int grid, double* x) const {
    if (grid == 0) {
        for (int i = 0; i < ncols_; i++) {
            x[i] = xll_ + i * cellsize_;
        }
    }
}

void BmiOpenLISEM::GetGridY(int grid, double* y) const {
    if (grid == 0) {
        for (int i = 0; i < nrows_; i++) {
            y[i] = yll_ + i * cellsize_;
        }
    }
}

void BmiOpenLISEM::InitializeVariableRegistry() {
    // Register key input/output variables
    // Format: name, type, units, location, grid_id, is_input, is_output

    // Precipitation (input)
    RegisterVariable("rainfall_intensity", "double", "mm/h", "node", 0, true, false);

    // Surface water (output)
    RegisterVariable("surface_water_depth", "double", "m", "node", 0, false, true);
    RegisterVariable("surface_runoff", "double", "m3/s", "node", 0, false, true);
    RegisterVariable("surface_velocity", "double", "m/s", "node", 0, false, true);

    // Infiltration (output)
    RegisterVariable("infiltration_rate", "double", "mm/h", "node", 0, false, true);
    RegisterVariable("cumulative_infiltration", "double", "mm", "node", 0, false, true);

    // Soil moisture (output)
    RegisterVariable("soil_moisture_content", "double", "m3/m3", "node", 0, false, true);

    // Erosion (output)
    RegisterVariable("detachment_rate", "double", "kg/m2/s", "node", 0, false, true);
    RegisterVariable("sediment_concentration", "double", "kg/m3", "node", 0, false, true);
    RegisterVariable("cumulative_erosion", "double", "kg/m2", "node", 0, false, true);

    // Channel flow (output)
    RegisterVariable("channel_discharge", "double", "m3/s", "node", 0, false, true);
    RegisterVariable("channel_water_depth", "double", "m", "node", 0, false, true);

    // Update variable pointers
    UpdateVariablePointers();
}

void BmiOpenLISEM::RegisterVariable(const std::string& name,
                                     const std::string& type,
                                     const std::string& units,
                                     const std::string& location,
                                     int grid_id,
                                     bool is_input,
                                     bool is_output) {
    VarInfo info;
    info.type = type;
    info.units = units;
    info.location = location;
    info.grid_id = grid_id;
    info.data_ptr = nullptr;

    var_info_[name] = info;

    if (is_input) {
        input_var_names_.push_back(name);
    }
    if (is_output) {
        output_var_names_.push_back(name);
    }
}

void* BmiOpenLISEM::GetVariablePointer(const std::string& name) {
    if (!model_) return nullptr;

    auto it = var_info_.find(name);
    if (it == var_info_.end()) return nullptr;

    // Return pointer to the actual model data
    // Map variable names to TWorld class members
    if (name == "surface_water_depth" && model_->WH) {
        return model_->WH->data[0];
    }
    else if (name == "surface_runoff" && model_->Q) {
        return model_->Q->data[0];
    }
    else if (name == "surface_velocity" && model_->V) {
        return model_->V->data[0];
    }
    else if (name == "infiltration_rate" && model_->Infil) {
        return model_->Infil->data[0];
    }
    else if (name == "cumulative_infiltration" && model_->InfilCum) {
        return model_->InfilCum->data[0];
    }
    else if (name == "soil_moisture_content" && model_->Theta1) {
        return model_->Theta1->data[0];
    }
    else if (name == "detachment_rate" && model_->DEP) {
        return model_->DEP->data[0];
    }
    else if (name == "sediment_concentration" && model_->Conc) {
        return model_->Conc->data[0];
    }
    else if (name == "cumulative_erosion" && model_->SedCum) {
        return model_->SedCum->data[0];
    }
    else if (name == "channel_discharge" && model_->ChannelQ) {
        return model_->ChannelQ->data[0];
    }
    else if (name == "channel_water_depth" && model_->ChannelWH) {
        return model_->ChannelWH->data[0];
    }
    else if (name == "rainfall_intensity" && model_->Rainc) {
        return model_->Rainc->data[0];
    }

    return nullptr;
}

void BmiOpenLISEM::UpdateVariablePointers() {
    // Update all registered variable pointers
    for (auto& pair : var_info_) {
        pair.second.data_ptr = GetVariablePointer(pair.first);
    }
}

} // namespace bmi

/* C interface implementation */

extern "C" {

void* bmi_new() {
    return new bmi::BmiOpenLISEM();
}

int bmi_initialize(void* handle, const char* config_file) {
    if (!handle) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    return model->Initialize(std::string(config_file));
}

int bmi_update(void* handle) {
    if (!handle) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    return model->Update();
}

int bmi_update_until(void* handle, double time) {
    if (!handle) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    return model->UpdateUntil(time);
}

int bmi_finalize(void* handle) {
    if (!handle) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    return model->Finalize();
}

int bmi_get_current_time(void* handle, double* time) {
    if (!handle || !time) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    *time = model->GetCurrentTime();
    return BMI_SUCCESS;
}

int bmi_get_start_time(void* handle, double* time) {
    if (!handle || !time) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    *time = model->GetStartTime();
    return BMI_SUCCESS;
}

int bmi_get_end_time(void* handle, double* time) {
    if (!handle || !time) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    *time = model->GetEndTime();
    return BMI_SUCCESS;
}

int bmi_get_time_step(void* handle, double* dt) {
    if (!handle || !dt) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    *dt = model->GetTimeStep();
    return BMI_SUCCESS;
}

int bmi_get_value(void* handle, const char* name, void* dest) {
    if (!handle || !name || !dest) return BMI_FAILURE;
    try {
        bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
        model->GetValue(std::string(name), dest);
        return BMI_SUCCESS;
    }
    catch (...) {
        return BMI_FAILURE;
    }
}

int bmi_set_value(void* handle, const char* name, void* src) {
    if (!handle || !name || !src) return BMI_FAILURE;
    try {
        bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
        model->SetValue(std::string(name), src);
        return BMI_SUCCESS;
    }
    catch (...) {
        return BMI_FAILURE;
    }
}

int bmi_get_grid_size(void* handle, int grid, int* size) {
    if (!handle || !size) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    *size = model->GetGridSize(grid);
    return BMI_SUCCESS;
}

int bmi_get_grid_rank(void* handle, int grid, int* rank) {
    if (!handle || !rank) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    *rank = model->GetGridRank(grid);
    return BMI_SUCCESS;
}

int bmi_get_grid_shape(void* handle, int grid, int* shape) {
    if (!handle || !shape) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    model->GetGridShape(grid, shape);
    return BMI_SUCCESS;
}

int bmi_get_grid_spacing(void* handle, int grid, double* spacing) {
    if (!handle || !spacing) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    model->GetGridSpacing(grid, spacing);
    return BMI_SUCCESS;
}

int bmi_get_grid_origin(void* handle, int grid, double* origin) {
    if (!handle || !origin) return BMI_FAILURE;
    bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
    model->GetGridOrigin(grid, origin);
    return BMI_SUCCESS;
}

void bmi_delete(void* handle) {
    if (handle) {
        bmi::BmiOpenLISEM* model = static_cast<bmi::BmiOpenLISEM*>(handle);
        delete model;
    }
}

} // extern "C"
