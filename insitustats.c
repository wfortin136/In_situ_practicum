#include "insitustats.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <mpi.h>
#include <stdint.h>

int* set_indicies(int i_x, int i_y, int i_z)
{
  int* x = (int*)malloc(sizeof(int)*3);
  if(!x){
    perror("malloc");
  }

  x[0]=i_x;
  x[1]=i_y;
  x[2]=i_z;

  return x;
}

field_val* field_val_new(char* name, double* data_set){
  
  field_val* f = (field_val*)malloc(sizeof(field_val));
  if(!f){
    perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }
/*
  char* n = (char*)malloc(strlen(name)*sizeof(char));
  if(!n){
    perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }
  */
  f->name = name;
  f->field_str = str_stat_new(data_set);

  return f;

}

var_stats* str_stat_new(double* data_set)
{

  var_stats * n = (var_stats*)malloc(sizeof(var_stats));
  if(!n){
    perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }

  n->histo = str_histo_new();
  n->data_set = data_set;
  n->mean = 0;
  n->min = 0;
  n->max =0;
  
  return n;
}

//Different then str_stat_new because it doesn't pass data set
//May consider passing empty data set instead of allocating
//within function
var_stats* str_stat_new_g()
{

  var_stats * n = (var_stats*)malloc(sizeof(var_stats));
  if(!n){
    perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }

  n->histo = str_histo_new();
  n->mean = 0;
  n->min = 0;
  n->max =0;
  
  return n;
}

histogram* str_histo_new()
{
  histogram* h = (histogram*)malloc(sizeof(histogram));
  if(!(h)){
    perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }

  h->count = (int*)malloc(sizeof(int));
  if(!(h->count)){
    perror("malloc");
  }

  h->index = (int*)malloc(sizeof(int));
  if(!(h->index)){
    perror("malloc");
  }

  return h;
}

void compute_stats(var_stats * var_struct,int* ind)
{
  unsigned long long i, index;
  double sum=0, min, max;

  index = ind[0] * ind[1] * ind[2];//l_x, l_y, and l_z
  min= var_struct->data_set[0];
  max= var_struct->data_set[0];
  
  for(i=0; i<index; i++){
    //printf("j:%llu, k:%llu, i:%llu, index:%llu\n",j,k,i,index);

    sum+= var_struct->data_set[i];  
    if(var_struct->data_set[i]>max) max=var_struct->data_set[i];
    if(var_struct->data_set[i]<min) min=var_struct->data_set[i];
 
    //printf("index:%llu, volume:%f \n",i, var_struct->data_set[i]);
  }
  var_struct->mean= (sum/index);
  var_struct->min=min;
  var_struct->max=max;
  //Pulling out compute_histo into made program because we will compute
  //using the global min and max

  //compute_histo(var_struct->histo, var_struct->data_set, min, max); 
}


//static void compute_histo(var_stats * var_struct)
void compute_histo(histogram * histo, double* data_set, 
    double min, double max, int* ind)
{
  int i, index, lower, upper;
  int j, range, step, numcols;

  index = ind[0] * ind[1] * ind[2];//l_x, l_y, and l_z

  //Try and determine reasonable steps and ranges for histogram
  //make sure max int is larger than double max
  range = (max)-(min)+1;
  
  if((2 % range)!=0) range+=1;//make sure range is even
  numcols=range;
  step=range/numcols;  
  //For now, assume more than 10 columns is too much. 
  //Arbitrary at this point
  while((numcols>10)){
    numcols/=2;
    //if((2 % range)==1) range+=1;//turn number even if odd
    step=range/numcols;//take integer portion of division

    //printf("range: %i, step:%i,numcols:%i \n", range, step, numcols);
  }
  //printf("range: %i, step:%i,numcols:%i \n", range, step, numcols);
  histo->numcols=numcols;
  
  //reseize array
  histo->count = realloc(histo->count,sizeof(int)*numcols);
  if(!(histo->count)){
      perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }
  
  //resize array
  histo->index = realloc(histo->index,sizeof(int)*(numcols+1));
  if(!(histo->index)){
      perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }
  histo->index[0]=min;

  for(j=0; j<numcols; j++){
    histo->index[j+1]=histo->index[j]+step;
    histo->count[j]=0;
  }

  //parse each value within data set and place in correct bin
  for(i=0; i<index; i++){
    
    for(j=0; j<numcols; j++){
      //printf("j:%d, i:%d, index:%d\n",j,i,index);
      lower= histo->index[j];
      upper= histo->index[j+1];
      //printf("L:%d U:%d j:%d\n", lower,upper, j);
      //printf("data: %f\n",data_set[i]);
      if((data_set[i]>=lower) && (data_set[i]<upper)){
        histo->count[j]+=1;
        //printf("data: %d\n",j);
      }

      //printf("data: %d\n",j);
    }
    //printf("index:%llu, volume:%f \n",i, var_struct->data_set[i]);
  }

}

void print_stats(char* var_name, var_stats * var_struct)
{
  int i;

  printf("****%s**** \nMean: %f\nMin: %f\nMax: %f\n", var_name, var_struct->mean,
				  var_struct->min, var_struct->max);  
  for(i=0; i<strlen(var_name)+8; i++)
  {
    printf("*");
  }
  printf("\n");
}

void print_histo(char* var_name, histogram* histo)
{
  int i, x;
  int lower, upper, count;
  
  x=histo->numcols;
  printf("***Histo-%s***\n", var_name);
  for(i=0; i<x; i++){
    upper=histo->index[i+1];
    lower=histo->index[i];
    count=histo->count[i]; 
    printf("%i-%i: %i\n", lower,upper, count);
  }
}

void free_indicies(int* i)
{
  free(i);
}

void free_field_val(field_val* fv)
{
  free_var_struct(fv->field_str);
  //free(fv->name);
  free(fv);
}

void free_var_struct(var_stats * var_struct)
{
  free_histo_struct(var_struct->histo);
  free(var_struct);
}

void free_histo_struct(histogram * histo)
{
  free(histo->index);
  free(histo->count);
  free(histo);
}
