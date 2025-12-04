/*************************************************************************
**  BMI (Basic Model Interface) header for openLISEM
**  Based on CSDMS BMI 2.0 specification
**  https://bmi.readthedocs.io/
**************************************************************************/

#ifndef BMI_H
#define BMI_H

#ifdef __cplusplus
extern "C" {
#endif

/* BMI return codes */
#define BMI_SUCCESS (0)
#define BMI_FAILURE (1)

/* BMI function declarations */

/* Model control functions */
int initialize(const char *config_file);
int update();
int update_until(double time);
int finalize();

/* Model information functions */
void get_component_name(char *name);
int get_input_item_count();
int get_output_item_count();
void get_input_var_names(char **names);
void get_output_var_names(char **names);

/* Variable information functions */
int get_var_grid(const char *name);
void get_var_type(const char *name, char *type);
void get_var_units(const char *name, char *units);
int get_var_itemsize(const char *name);
int get_var_nbytes(const char *name);
void get_var_location(const char *name, char *location);

/* Time functions */
double get_current_time();
double get_start_time();
double get_end_time();
double get_time_step();
void get_time_units(char *units);

/* Variable getter and setter functions */
int get_value(const char *name, void *dest);
int get_value_ptr(const char *name, void **dest_ptr);
int get_value_at_indices(const char *name, void *dest, int *inds, int count);
int set_value(const char *name, void *src);
int set_value_at_indices(const char *name, int *inds, int count, void *src);

/* Grid information functions */
int get_grid_rank(int grid);
int get_grid_size(int grid);
void get_grid_type(int grid, char *type);
int get_grid_shape(int grid, int *shape);
int get_grid_spacing(int grid, double *spacing);
int get_grid_origin(int grid, double *origin);
int get_grid_x(int grid, double *x);
int get_grid_y(int grid, double *y);
int get_grid_z(int grid, double *z);
int get_grid_node_count(int grid);
int get_grid_edge_count(int grid);
int get_grid_face_count(int grid);
int get_grid_edge_nodes(int grid, int *edge_nodes);
int get_grid_face_edges(int grid, int *face_edges);
int get_grid_face_nodes(int grid, int *face_nodes);
int get_grid_nodes_per_face(int grid, int *nodes_per_face);

#ifdef __cplusplus
}
#endif

#endif /* BMI_H */
