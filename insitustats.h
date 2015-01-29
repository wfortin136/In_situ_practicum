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
  //int* histo_index;//index values for ranges in histogram
  //int* histo;//count for index values
  //int numcols; //number of columns in histogram
};
typedef struct var_stats var_stats;
//********************

var_stats* str_stat_new(double* volume);

var_stats* str_stat_new_g();

histogram* str_histo_new();

int* set_indicies(int i_x, int i_y, int i_z);

void free_indicies(int* i);

void compute_stats(var_stats * var_struct, int* ind);

void compute_histo(histogram * histo, double* data_set, 
	double min, double max, int* ind);

void print_stats(char* var_name, var_stats * var_struct);

void print_histo(char* var_name, histogram* histo);

void free_var_struct(var_stats * var_struct);

void free_histo_struct(histogram * histo);
