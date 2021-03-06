#pragma once

//****************
//Struct for statistics
struct histogram
{
  int* index;//index values for ranges in histogram
  int* count;//count for index values
  int numcols; //number of columns in histogram
};
typedef struct histogram histogram;

struct var_stats
{
  double* data_set;
  double mean; 
  double min;
  double max;
  histogram* histo;
};
typedef struct var_stats var_stats;

struct field_val
{
  char* name;
  var_stats* field_str;
};
typedef struct field_val field_val;

/*
struct field_vals_array
{
  field_val* values;
  int array_size;
};
typedef struct field_vals_array field_vals_array;
*/

//typedef struct field_val* two_d_slices;
typedef struct field_val** vals_array;

struct two_d_slices
{
  field_val** slices;
  int num_slices;
};
typedef struct two_d_slices two_d_slices;

//********************

vals_array new_field_vals_array();

vals_array new_val_in_field_array(vals_array array, char* name, double* data_set);

field_val* field_val_new(char* name, double* data_set);

field_val* field_val_new_empty(char* name);

two_d_slices* create_2d_slices(char* name, double* data_set, int* indicies, int dim);

histogram* str_histo_new();

int* set_indicies(int i_x, int i_y, int i_z);

void free_indicies(int* i);

void compute_var_stats(field_val* f_val, int* indicies);

char* get_var_name(field_val* f_val);

double* get_mean_ptr(field_val* f_val);

double* get_min_ptr(field_val* f_val);

double* get_max_ptr(field_val* f_val);

void calc_histogram(field_val* fv_histo, field_val* fv_data, 
    double min, double max, int* indicies);

int* get_histo_array(field_val* fv);

int get_histo_numcols(field_val* fv);

void print_stats(char* var_name, var_stats * var_struct);

void print_histo(char* var_name, histogram* histo);

void free_field_val(field_val* fv);

void free_2d_slices(two_d_slices* slices_str);
