/*************************************************************************
**  Python bindings for BMI openLISEM using pybind11
**************************************************************************/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "bmi_openlisem.h"

namespace py = pybind11;

PYBIND11_MODULE(bmi_openlisem, m) {
    m.doc() = "BMI wrapper for openLISEM hydrological and erosion model";

    py::class_<bmi::BmiOpenLISEM>(m, "BmiOpenLISEM")
        .def(py::init<>())

        // Control functions
        .def("initialize", &bmi::BmiOpenLISEM::Initialize,
             py::arg("config_file"),
             "Initialize the model with a configuration file (runfile)")

        .def("update", &bmi::BmiOpenLISEM::Update,
             "Advance the model by one time step")

        .def("update_until", &bmi::BmiOpenLISEM::UpdateUntil,
             py::arg("time"),
             "Advance the model until the given time")

        .def("finalize", &bmi::BmiOpenLISEM::Finalize,
             "Finalize the model")

        // Model information
        .def("get_component_name", &bmi::BmiOpenLISEM::GetComponentName,
             "Get the name of the model component")

        .def("get_input_item_count", &bmi::BmiOpenLISEM::GetInputItemCount,
             "Get the number of input variables")

        .def("get_output_item_count", &bmi::BmiOpenLISEM::GetOutputItemCount,
             "Get the number of output variables")

        .def("get_input_var_names", &bmi::BmiOpenLISEM::GetInputVarNames,
             "Get the list of input variable names")

        .def("get_output_var_names", &bmi::BmiOpenLISEM::GetOutputVarNames,
             "Get the list of output variable names")

        // Variable information
        .def("get_var_grid", &bmi::BmiOpenLISEM::GetVarGrid,
             py::arg("name"),
             "Get the grid identifier for a variable")

        .def("get_var_type", &bmi::BmiOpenLISEM::GetVarType,
             py::arg("name"),
             "Get the data type of a variable")

        .def("get_var_units", &bmi::BmiOpenLISEM::GetVarUnits,
             py::arg("name"),
             "Get the units of a variable")

        .def("get_var_itemsize", &bmi::BmiOpenLISEM::GetVarItemsize,
             py::arg("name"),
             "Get the size of a single element of a variable in bytes")

        .def("get_var_nbytes", &bmi::BmiOpenLISEM::GetVarNbytes,
             py::arg("name"),
             "Get the total size of a variable in bytes")

        .def("get_var_location", &bmi::BmiOpenLISEM::GetVarLocation,
             py::arg("name"),
             "Get the location (node, edge, face) of a variable")

        // Time functions
        .def("get_current_time", &bmi::BmiOpenLISEM::GetCurrentTime,
             "Get the current model time")

        .def("get_start_time", &bmi::BmiOpenLISEM::GetStartTime,
             "Get the start time of the simulation")

        .def("get_end_time", &bmi::BmiOpenLISEM::GetEndTime,
             "Get the end time of the simulation")

        .def("get_time_step", &bmi::BmiOpenLISEM::GetTimeStep,
             "Get the model time step")

        .def("get_time_units", &bmi::BmiOpenLISEM::GetTimeUnits,
             "Get the time units")

        // Variable getters and setters with numpy array support
        .def("get_value", [](bmi::BmiOpenLISEM& self, const std::string& name) {
            int grid = self.GetVarGrid(name);
            if (grid == 0) {
                int shape[2];
                self.GetGridShape(grid, shape);

                auto result = py::array_t<double>({shape[0], shape[1]});
                self.GetValue(name, result.mutable_data());
                return result;
            }
            throw std::runtime_error("Variable not found or invalid grid");
        }, py::arg("name"), "Get the value of a variable as a numpy array")

        .def("set_value", [](bmi::BmiOpenLISEM& self, const std::string& name,
                            py::array_t<double> values) {
            self.SetValue(name, values.mutable_data());
        }, py::arg("name"), py::arg("values"),
           "Set the value of a variable from a numpy array")

        .def("get_value_ptr", [](bmi::BmiOpenLISEM& self, const std::string& name) {
            void* ptr = self.GetValuePtr(name);
            if (!ptr) {
                throw std::runtime_error("Variable not found: " + name);
            }

            int grid = self.GetVarGrid(name);
            if (grid == 0) {
                int shape[2];
                self.GetGridShape(grid, shape);

                // Return a numpy array view (no copy)
                return py::array_t<double>(
                    {shape[0], shape[1]},
                    {shape[1] * sizeof(double), sizeof(double)},
                    static_cast<double*>(ptr),
                    py::cast(self)
                );
            }
            throw std::runtime_error("Invalid grid");
        }, py::arg("name"), "Get a direct pointer to a variable as a numpy array view")

        // Grid information
        .def("get_grid_rank", &bmi::BmiOpenLISEM::GetGridRank,
             py::arg("grid"),
             "Get the rank (number of dimensions) of a grid")

        .def("get_grid_size", &bmi::BmiOpenLISEM::GetGridSize,
             py::arg("grid"),
             "Get the total number of elements in a grid")

        .def("get_grid_type", &bmi::BmiOpenLISEM::GetGridType,
             py::arg("grid"),
             "Get the type of a grid")

        .def("get_grid_shape", [](bmi::BmiOpenLISEM& self, int grid) {
            int shape[2];
            self.GetGridShape(grid, shape);
            return py::make_tuple(shape[0], shape[1]);
        }, py::arg("grid"), "Get the shape of a grid")

        .def("get_grid_spacing", [](bmi::BmiOpenLISEM& self, int grid) {
            double spacing[2];
            self.GetGridSpacing(grid, spacing);
            return py::make_tuple(spacing[0], spacing[1]);
        }, py::arg("grid"), "Get the spacing of a grid")

        .def("get_grid_origin", [](bmi::BmiOpenLISEM& self, int grid) {
            double origin[2];
            self.GetGridOrigin(grid, origin);
            return py::make_tuple(origin[0], origin[1]);
        }, py::arg("grid"), "Get the origin of a grid");

    // Module-level constants
    m.attr("BMI_SUCCESS") = py::int_(BMI_SUCCESS);
    m.attr("BMI_FAILURE") = py::int_(BMI_FAILURE);
}
