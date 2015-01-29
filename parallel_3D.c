//
// Parallel_3D_Volume.c
//  
//
//  Created by Venkatram Vishwanath on 12/17/14.
//  Updated by Billy Fortin on 1/23/2015
//

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <mpi.h>
//Billy added
#include <stdint.h>

int dim = 3;

// global dimensions of 2D volume
static int g_x = 0;
static int g_y = 0;
static int g_z = 0;

// per-process dimensions of each sub-block
static int l_x = 0;
static int l_y = 0;
static int l_z = 0;

//****************
//Billy
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

//

static int parse_args(int argc, char **argv);
static void usage(void);

// Assigns value based on global index
static void gen_data_global_index (double* volume,
                     long long my_off_z, long long my_off_y, long long my_off_x);

static void gen_data (double* volume,
                     long long my_off_z, long long my_off_y, long long my_off_x);

static void gen_data_sparse (double* volume,
                     long long my_off_z, long long my_off_y, long long my_off_x);

//**************************************
//Billy
var_stats* str_stat_new(double* volume);

var_stats* str_stat_new_g();

histogram* str_histo_new();

static void compute_stats(var_stats * var_struct);

//static void compute_histo(var_stats * var_struct);
static void compute_histo(histogram * histo, double* data_set, double min, double max);

static void print_stats(char* var_name, var_stats * var_struct);

static void print_histo(char* var_name, histogram* histo);

static void free_var_struct(var_stats * var_struct);

static void free_histo_struct(histogram * histo);
//****************************************

//

int main(int argc, char **argv)
{
  int nprocs, rank, rc;
  int ret;
  uint64_t  index, slice;
  int tot_blocks_x, tot_blocks_y, tot_blocks_z, tot_blocks;
  uint64_t start_extents_z, start_extents_y, start_extents_x;
  int block_id_x, block_id_y, block_id_z;
  int i,j,k;
  
 
  // The buffers/ variables
  // Let's have Pressure, Temperature, Density
  // TODO: Make the num of variables a command line argument
  double* pressure = 0;
  double* temperature = 0;
  double* density = 0;

  //***************************************************
  //Billy
  double mean_2=0;
  double mean_1=0;
  double mean =0; //was using mean to test what was causing seg fault
  //mean = (double*) malloc(sizeof(double));
  //****************************************************

  //mean_1 = (double*) malloc(sizeof(double));
  //mean_2 = (double*) malloc(sizeof(double));
  
  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //sleep(60);
  ////////////////////////////
  // parse args on rank 0
  ////////////////////////////
  if(rank == 0){
    
      ret = parse_args(argc, argv);
      if(ret < 0){
          usage();
          MPI_Abort(MPI_COMM_WORLD, -1);
      }
      
      // check if the num procs is appropriate
      tot_blocks = (g_z/l_z) * (g_y/l_y) * (g_x/l_x);
    
      if(tot_blocks != nprocs){
          printf("Error: number of blocks (%d) doesn't match\
				number of procs (%d)\n", tot_blocks, nprocs);
          MPI_Abort(MPI_COMM_WORLD, -1);
      }

      //***************************************************
      //Billy
      //mean_2 = (double*) malloc(nprocs*sizeof(double));
      //***********************************************
  }
  
  /////////////////////////////
  // share the command line args and other params
  /////////////////////////////
  
  MPI_Bcast(&g_z, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&g_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&l_z, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&l_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&l_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&tot_blocks, 1, MPI_INT, 0, MPI_COMM_WORLD);

 
  // figure out some kind of block distribution
  tot_blocks_z = (g_z/l_z);
  tot_blocks_y = (g_y/l_y);
  tot_blocks_x = (g_x/l_x);
  
  // start extents in Z, Y, X for my block
  if (nprocs == 1)
  {
    start_extents_z = 0; start_extents_y = 0; start_extents_x = 0;
  }
  else
  {
    start_extents_z = (rank / (tot_blocks_y * tot_blocks_x)) * l_z;
    slice = rank % (tot_blocks_y * tot_blocks_x);
    start_extents_y = (slice / tot_blocks_x) * l_y;
    start_extents_x = (slice % tot_blocks_x) * l_x;
  }
  
  
  block_id_z = start_extents_z / l_z;
  block_id_y = start_extents_y / l_y;
  block_id_x = start_extents_x / l_x;
 
  // Print Info
  if (0 == rank){
    printf("Global Dimensions %dX%dX%d: Local Dimensions %dX%dX%d \n", \
          g_z, g_y, g_x, l_z, l_y, l_x);
    printf("Total Blocks are %dX%dX%d \n", tot_blocks_z, tot_blocks_y, tot_blocks_x);
  }

  //////////////////////////////////
  // allocate the variables and
  // intialize to a pattern
  //////////////////////////////////
  
  // Variable allocation
  pressure = (double*) malloc (sizeof(double) * l_z * l_y *l_x);
  if(!pressure){
      perror("malloc");
      MPI_Abort(MPI_COMM_WORLD, -1);
  }
  
  temperature = (double*) malloc (sizeof(double) * l_z * l_y *l_x);
  if(!temperature){
    perror("malloc");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  
  density = (double*) malloc (sizeof(double) * l_z * l_y *l_x);
  if(!density){
    perror("malloc");
    MPI_Abort(MPI_COMM_WORLD, -1);
  }
  
  
  // Initialize Variables to a value for testing
  // Note: Currently Set this to the Global Index Value
 /*
  for(k=0; k<l_z; k++){
    for(j=0; j<l_y; j++){
      for(i=0; i<l_x; i++){
        index = (l_x * l_y * k) + (l_x*j) + i;
        pressure[index] = (g_x * g_y * (start_extents_z + k))
                       + (g_x * (start_extents_y + j)) + start_extents_x + i;
        
        temperature[index] = (g_x * g_y * (start_extents_z + k))
                          + (g_x * (start_extents_y + j)) + start_extents_x + i;
       
        density[index] = (g_x * g_y * (start_extents_z + k))
                        + (g_x * (start_extents_y + j)) + start_extents_x + i;

      }
    }
  }
 */ 

  // Intialize the variables to the values of the global index
  gen_data_global_index (pressure, start_extents_z, start_extents_y, start_extents_x);
  gen_data_global_index (temperature, start_extents_z, start_extents_y, start_extents_x);
  gen_data_global_index (density, start_extents_z, start_extents_y, start_extents_x);

  // DEBUG: Print the values of the variables..
    
  //MPI_Barrier(MPI_COMM_WORLD);
  
  /////////////////////////////
  // Iterate over multiple timesteps?
  // Compute several analyses?
  /////////////////////////////
  
  //**********************************  
  //Billy
  //Define struct for stats
  var_stats* l_pressure_str = str_stat_new(pressure);
  var_stats* l_temp_str = str_stat_new(temperature);
  var_stats* l_density_str = str_stat_new(density);
  
  var_stats* g_pressure_str = str_stat_new_g();
  var_stats* g_temp_str = str_stat_new_g();
  var_stats* g_density_str = str_stat_new_g();
/*
  for(j=0; j<nprocs; j++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==j){
      //MPI_Barrier(MPI_COMM_WORLD); 
      for(i=0; i<(l_x * l_y *l_z); i++)
      {
	printf("Index:%i, Volume:%f, Dataset:%f\n",i, pressure[i],pressure_str->data_set[i]);  
      }
    }
  }
*/
  compute_stats(l_pressure_str);
  compute_stats(l_temp_str);
  compute_stats(l_density_str);
  mean =(double)rank;

  /*
  print_histo("Pressure", pressure_str->histo);
  print_histo("Temperature", temp_str->histo);
  print_histo("Density", density_str->histo);
  print_stats("Pressure", pressure_str);
  print_stats("Temperature", temp_str);
  print_stats("Density", density_str);
  */
  
  //
  //*************************************************
  //MPI_Barrier(MPI_COMM_WORLD);
 
  ////////////////////////////////
  //Collate Statistics across Processes
  ////////////////////////////////
  //************************************************
  int s=0;
 
  //MPI_Type_size(MPI_DOUBLE, &s);
  //printf("System: %li MPI: %d\n", sizeof(double), s);
  //Billy    
  if(rank==0){
    printf("**BEFORE**\n");
    printf("Pressure\n...Glebal Mean:%f\n...Global Min:%f\n...Global Max:%f\n", 
        g_pressure_str->mean, g_pressure_str->min, g_pressure_str->max);
    printf("Density\n...Glebal Mean:%f\n...Global Min:%f\n...Global Max:%f\n", 
        g_density_str->mean, g_density_str->min, g_density_str->max);
    printf("Temperature\n...Glebal Mean:%f\n...Global Min:%f\n...Global Max:%f\n", 
        g_temp_str->mean, g_temp_str->min, g_temp_str->max);
  }

  //printf("mean: %f \n", pressure_str->mean);

//**************TEST WITH NOT STRUCT LOCAL MEAN****************
//
  //MPI_Barrier(MPI_COMM_WORLD);
  //MPI_Allreduce(&mean_1, &mean_2, 1,MPI_DOUBLE, 
  //    MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&(l_pressure_str->mean), &(g_pressure_str->mean), 1,MPI_DOUBLE, 
      MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&(l_density_str->mean), &(g_density_str->mean), 1,MPI_DOUBLE, 
      MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&(l_temp_str->mean), &(g_temp_str->mean), 1,MPI_DOUBLE, 
      MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce(&(l_pressure_str->min), &(g_pressure_str->min), 1,MPI_DOUBLE, 
      MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&(l_density_str->min), &(g_density_str->min), 1,MPI_DOUBLE, 
      MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce(&(l_temp_str->min), &(g_temp_str->min), 1,MPI_DOUBLE, 
      MPI_MIN, MPI_COMM_WORLD);
  
  MPI_Allreduce(&(l_pressure_str->max), &(g_pressure_str->max), 1,MPI_DOUBLE, 
      MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&(l_density_str->max), &(g_density_str->max), 1,MPI_DOUBLE, 
      MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&(l_temp_str->max), &(g_temp_str->max), 1,MPI_DOUBLE, 
      MPI_MAX, MPI_COMM_WORLD);
 
  g_pressure_str->mean/=nprocs;
  g_density_str->mean/=nprocs;
  g_temp_str->mean/=nprocs;

 
  compute_histo(l_pressure_str->histo, l_pressure_str->data_set, 
      g_pressure_str->min, g_pressure_str->max);
  compute_histo(l_density_str->histo, l_density_str->data_set, 
      g_density_str->min, g_density_str->max);
  compute_histo(l_temp_str->histo, l_temp_str->data_set, 
      g_temp_str->min, g_temp_str->max); 
  
  //Need to compute global histogram with local data in order to properly
  //allocate space for MPI_reduce
  //if(rank==0){
    compute_histo(g_pressure_str->histo, l_pressure_str->data_set, 
        g_pressure_str->min, g_pressure_str->max);
    compute_histo(g_density_str->histo, l_density_str->data_set, 
        g_density_str->min, g_density_str->max);
    compute_histo(g_temp_str->histo, l_temp_str->data_set, 
        g_temp_str->min, g_temp_str->max);
  //}
   
  MPI_Allreduce(l_pressure_str->histo->count, g_pressure_str->histo->count, 
      g_pressure_str->histo->numcols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l_temp_str->histo->count, g_temp_str->histo->count, 
      g_temp_str->histo->numcols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l_density_str->histo->count, g_density_str->histo->count, 
      g_density_str->histo->numcols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  
  //pressure_str->histo=

/*  
  for(i=0; i<; i++){
    
    compute_histo( 
  }
*/
  

  if(rank==0){
    //print_histo("Pressure", l_pressure_str->histo);
    printf("**AFTER**\n");
    printf("Pressure\n...Glebal Mean:%f\n...Global Min:%f\n...Global Max:%f\n", 
        g_pressure_str->mean, g_pressure_str->min, g_pressure_str->max);
    printf("Density\n...Glebal Mean:%f\n...Global Min:%f\n...Global Max:%f\n", 
        g_density_str->mean, g_density_str->min, g_density_str->max);
    printf("Temperature\n...Glebal Mean:%f\n...Global Min:%f\n...Global Max:%f\n", 
        g_temp_str->mean, g_temp_str->min, g_temp_str->max);
    print_histo("Pressure", g_pressure_str->histo);
    print_histo("Temperature", g_temp_str->histo);
    print_histo("Density", g_density_str->histo);
    //printf("mean test: %f\n", *mean_2);
  }
  
  //
  //********************************************************
  
  //MPI_Barrier(MPI_COMM_WORLD);
  /////////////////////////////
  // Clean up heap variables
  /////////////////////////////

  //***************************************
  //Billy
  //free(mean);
  //free(mean_1); 
  //free(mean_2);
  
  free_var_struct(l_pressure_str);
  free_var_struct(l_temp_str);
  free_var_struct(l_density_str);
  free_var_struct(g_pressure_str);
  free_var_struct(g_temp_str);
  free_var_struct(g_density_str);
  //
  
  if (pressure){
    free(pressure);
    pressure = 0;
  }

  if (temperature){
    free(temperature);
    temperature = 0;
  }
  
  if (density){
    free(density);
    density = 0;
  }
  fflush(stdout);
  //
  //printf("test1\n");
  //MPI_Barrier(MPI_COMM_WORLD);
  //****************************************
  MPI_Finalize();
  //printf("test: %d\n", rc);
  return(0);
}


/* parse_args()
 *
 * parses command line argument
 *
 * returns 0 on success, -1 on failure
 */
static int parse_args(int argc, char **argv)
{
    char flags[] = "g:l:";
    int one_opt = 0;
    
  while((one_opt = getopt(argc, argv, flags)) != EOF){
  // postpone error checking for after while loop */
    switch(one_opt){
        case('g'):
            sscanf(optarg, "%dx%dx%d", &g_z, &g_y, &g_x);
            break;
        case('l'):
            sscanf(optarg, "%dx%dx%d", &l_z, &l_y, &l_x);
            break;
        case('?'):
            return(-1);
    }
  }
  
    //printf ("Global Values : %d %d %d \n", g_z, g_y, g_x );
    //printf ("Local Values : %d %d %d \n", l_z, l_y, l_x );
  
  
  // need positive dimensions
  if(g_z < 1 || g_y < 1 || g_x < 1 ||l_z < 1 || l_y < 1 || l_x < 1 ) {
      printf("Error: bad dimension specification.\n");
      return(-1);
  }
  
  // need everything to be divisible
  if((g_z % l_z) || (g_y % l_y) || (g_x % l_x)){
      printf("Error: global dimensions and local dimensions aren't evenly divisible\n");
      return(-1);
  }

  return 0;
}



// prints usage instructions
static void usage(void)
{
  printf("Usage: <exec> 4x4x4 -l 2x2x2 \n");
  printf("  -g global dimensions\n");
  printf("  -l local (per-process) dimensions\n");
  printf("\n");
  printf("  test-one-side generates a 3D volume in parallel\n");
  
  return;
}

static void gen_data_global_index (double* volume,
                     long long my_off_z, long long my_off_y, long long my_off_x)
{
  unsigned long long i,j,k, index;

  for(k=0; k<l_z; k++){
    for(j=0; j<l_y; j++){
      for(i=0; i<l_x; i++){
        index = (l_x * l_y * k) + (l_x*j) + i;
        
	//printf("j:%llu, k:%llu, i:%llu, index:%llu\n",j,k,i,index);
 
        volume[index] = (g_x * g_y * (my_off_z + k)) +  \
                        (g_x * (my_off_y + j)) + my_off_x + i;
        //printf("index:%llu, volume:%f \n",index, volume[index]);
      }
    }
  }
  
}



/* generate a data set */
static void gen_data(double* volume,
                     long long my_off_z, long long my_off_y, long long my_off_x)
{
  unsigned long long i,j,k;
  double center[3];
  center[0] = g_x / 2.;
  center[1] = g_y / 2.;
  center[2] = g_z / 2.;
  
    //printf("l_x: %d, l_y: %d, l_z: %d\n", l_x, l_y, l_z);
  for(i = 0; i < l_z; i++)
  {
    double zdist = sin((i + my_off_z - center[2])/5.0)*center[2];
    //      float zdist = sinf((i + my_off_z - center[2])/g_y);
    for(j = 0; j < l_y; j++)
    {
      double ydist = sin((j + my_off_y - center[1])/3.0)*center[1];
        //      float ydist = sinf((j + my_off_y - center[1])/g_x);
      for(k = 0; k < l_x; k++)
      {
        double xdist = sin((k + my_off_x - center[0])/2.0)*center[0];
          //            float xdist = sinf((k + my_off_x - center[0])/g_z);
        volume[i * l_x * l_y + j * l_x + k] = sqrt(xdist * xdist + ydist *ydist + zdist * zdist);
        
      }
    }
  }
  
}

static void gen_data_sparse(double* volume,
                            long long my_off_z, long long my_off_y, long long my_off_x)
{
  unsigned long long i,j,k;
  double center[3];
  center[0] = g_x / 2.;
  center[1] = g_y / 2.;
  center[2] = g_z / 2.;
  
    //printf("l_x: %d, l_y: %d, l_z: %d\n", l_x, l_y, l_z);
  for(i = 0; i < l_z; i++)
  {
    double zdist = i + my_off_z - center[2];
    for(j = 0; j < l_y; j++)
    {
      double ydist = j + my_off_y - center[1];
      for(k = 0; k < l_x; k++)
      {
        double xdist = k + my_off_x - center[0];
        
        volume[i * l_x * l_y + j * l_x + k] = sqrt(xdist * xdist + ydist *ydist + zdist * zdist);
        
      }
    }
  }
  
  
}


//****************************************************
//Billy
var_stats* str_stat_new(double* volume)
{

  var_stats * n = (var_stats*)malloc(sizeof(var_stats));
  if(!n){
    perror("malloc");
      //MPI_Abort(MPI_COMM_WORLD, -1); SHOULD THIS BE INCLUDED?
  }

  n->histo = str_histo_new();
  n->data_set = volume;
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

  //Placeholder in case if empty data set for global structure
  /*
  n->data_set = (double*)malloc(sizeof(int));
  if(!(n->data_set)){
    perror("malloc");
  }
  */
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

static void compute_stats(var_stats * var_struct)
{
  unsigned long long i, index;
  double sum=0, min, max;
  
  index = l_x * l_y * l_z;
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
static void compute_histo(histogram * histo, double* data_set, double min, double max)
{
  int i, index, lower, upper;
  int j, range, step, numcols;

  index = l_x * l_y * l_z; 

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

static void print_stats(char* var_name, var_stats * var_struct)
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

static void print_histo(char* var_name, histogram* histo)
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
static void free_var_struct(var_stats * var_struct)
{
  free_histo_struct(var_struct->histo);
  free(var_struct);
}

static void free_histo_struct(histogram * histo)
{
  free(histo->index);
  free(histo->count);
  free(histo);
}
//*********************************************************
