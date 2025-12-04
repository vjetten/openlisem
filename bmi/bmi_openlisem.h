/*************************************************************************
**  BMI wrapper for openLISEM
**  Copyright (C) 2024
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**************************************************************************/

#ifndef BMI_OPENLISEM_H
#define BMI_OPENLISEM_H

#include <string>
#include <vector>
#include <map>
#include <memory>
#include "bmi.h"

// Forward declaration
class TWorld;
class QCoreApplication;

namespace bmi {

/**
 * @brief BMI wrapper class for openLISEM model
 *
 * This class implements the Basic Model Interface (BMI) for the openLISEM
 * hydrological and erosion model, enabling it to be used in model coupling
 * frameworks.
 */
class BmiOpenLISEM {
public:
    BmiOpenLISEM();
    ~BmiOpenLISEM();

    /* Model control functions */
    int Initialize(const std::string& config_file);
    int Update();
    int UpdateUntil(double time);
    int Finalize();

    /* Model information functions */
    std::string GetComponentName() const;
    int GetInputItemCount() const;
    int GetOutputItemCount() const;
    std::vector<std::string> GetInputVarNames() const;
    std::vector<std::string> GetOutputVarNames() const;

    /* Variable information functions */
    int GetVarGrid(const std::string& name) const;
    std::string GetVarType(const std::string& name) const;
    std::string GetVarUnits(const std::string& name) const;
    int GetVarItemsize(const std::string& name) const;
    int GetVarNbytes(const std::string& name) const;
    std::string GetVarLocation(const std::string& name) const;

    /* Time functions */
    double GetCurrentTime() const;
    double GetStartTime() const;
    double GetEndTime() const;
    double GetTimeStep() const;
    std::string GetTimeUnits() const;

    /* Variable getter and setter functions */
    void GetValue(const std::string& name, void* dest);
    void* GetValuePtr(const std::string& name);
    void GetValueAtIndices(const std::string& name, void* dest, int* inds, int count);
    void SetValue(const std::string& name, void* src);
    void SetValueAtIndices(const std::string& name, int* inds, int count, void* src);

    /* Grid information functions */
    int GetGridRank(int grid) const;
    int GetGridSize(int grid) const;
    std::string GetGridType(int grid) const;
    void GetGridShape(int grid, int* shape) const;
    void GetGridSpacing(int grid, double* spacing) const;
    void GetGridOrigin(int grid, double* origin) const;
    void GetGridX(int grid, double* x) const;
    void GetGridY(int grid, double* y) const;

private:
    std::unique_ptr<TWorld> model_;
    std::unique_ptr<QCoreApplication> qapp_;

    double current_time_;
    double start_time_;
    double end_time_;
    double time_step_;

    std::string config_file_;
    bool initialized_;
    bool finalized_;

    int nrows_;
    int ncols_;
    double cellsize_;
    double xll_;
    double yll_;

    // Variable name registry
    std::vector<std::string> input_var_names_;
    std::vector<std::string> output_var_names_;

    // Variable metadata storage
    struct VarInfo {
        std::string type;
        std::string units;
        std::string location;
        int grid_id;
        void* data_ptr;
    };

    std::map<std::string, VarInfo> var_info_;

    // Helper methods
    void InitializeVariableRegistry();
    void RegisterVariable(const std::string& name, const std::string& type,
                         const std::string& units, const std::string& location,
                         int grid_id, bool is_input, bool is_output);
    void* GetVariablePointer(const std::string& name);
    void UpdateVariablePointers();
};

} // namespace bmi

/* C interface wrapper functions */
extern "C" {
    void* bmi_new();
    int bmi_initialize(void* handle, const char* config_file);
    int bmi_update(void* handle);
    int bmi_update_until(void* handle, double time);
    int bmi_finalize(void* handle);
    int bmi_get_current_time(void* handle, double* time);
    int bmi_get_start_time(void* handle, double* time);
    int bmi_get_end_time(void* handle, double* time);
    int bmi_get_time_step(void* handle, double* dt);
    int bmi_get_value(void* handle, const char* name, void* dest);
    int bmi_set_value(void* handle, const char* name, void* src);
    int bmi_get_grid_size(void* handle, int grid, int* size);
    int bmi_get_grid_rank(void* handle, int grid, int* rank);
    int bmi_get_grid_shape(void* handle, int grid, int* shape);
    int bmi_get_grid_spacing(void* handle, int grid, double* spacing);
    int bmi_get_grid_origin(void* handle, int grid, double* origin);
    void bmi_delete(void* handle);
}

#endif /* BMI_OPENLISEM_H */
