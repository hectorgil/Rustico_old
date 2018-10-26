void get_skycuts_write_density_randoms(char *filename, double parameter_value[], double alpha, char *name_den_out);

void get_skycuts_write_density_data(char *filename, double parameter_value[],char *name_den_out);

void get_periodic_data(char *filename_data, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[]);

void get_skycuts_data(char *filename, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[], char *type_normalization_mode);

void get_skycuts_randoms(char *filename_data, double pos_x[], double pos_y[], double pos_z[], double weight[], double parameter_value[], char *type_normalization_mode, char *type_normalization_mode2);

void z_to_r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval);

long int get_number_used_lines_data(char *filename, double parameter_value[]);

long int get_number_used_lines_randoms(char *filename, double parameter_value[]);
