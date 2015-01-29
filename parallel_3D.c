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
#include <inttypes.h>
#include "insitustats.h"

int dim = 3;

// global dimensions of 2D volume
static int g_x = 0;
static int g_y = 0;
static int g_z = 0;

// per-process dimensions of each sub-block
static int l_x = 0;
static int l_y = 0;
static int l_z = 0;

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
  int* l_indicies=0; 
 
  double t1, t2, t3, t4, t5, t6;

  // The buffers/ variables
  // Let's have Pressure, Temperature, Density
  // TODO: Make the num of variables a command line argument
  double* pressure = 0;
  double* temperature = 0;
  double* density = 0;

  //************************
  //Set ranges for two of 3 dimensions

 
  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //sleep(60);
  ////////////////////////////
  // parse args on rank 0
  ////////////////////////////
  if(rank == 0){
    
      printf("Time Resolution: %E sec\n", MPI_Wtick());
      //1e-9 seconds is a billionth of a sec or nanosecond

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


  int m;
  for(m=0; m<nprocs; m++)
  {
    if(rank==m){
      printf("Rank: %d \n", rank);
      printf("id_x: %d, id_y: %d, id_z: %d\n", block_id_x, block_id_y, block_id_z);
      printf("id_x: %" PRIu64 ", id_y: %" PRIu64 ", id_z: %" PRIu64 "\n", start_extents_x, 
          start_extents_y, start_extents_z);
      fflush(stdout);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

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
  
  l_indicies = set_indicies(l_x, l_y, l_z);
  
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
  field_val* l_p = field_val_new("pressure", pressure);

  //
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
  compute_stats(l_pressure_str, l_indicies);
  compute_stats(l_temp_str, l_indicies);
  compute_stats(l_density_str, l_indicies);

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

/*  
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
*/
  //printf("mean: %f \n", pressure_str->mean);
  //collate_processes(l_pressure_str, g_pressure_str);
  //collate_processes("max",l_pressure_str);
  //collate_processes("min",l_pressure_str);
  
  t1=MPI_Wtime();
  
  MPI_Allreduce(&(l_pressure_str->mean), &(g_pressure_str->mean), 1,MPI_DOUBLE, 
      MPI_SUM, MPI_COMM_WORLD);
  t2=MPI_Wtime();

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
  t3=MPI_Wtime();
  g_pressure_str->mean/=nprocs;
  g_density_str->mean/=nprocs;
  g_temp_str->mean/=nprocs;

 
  compute_histo(l_pressure_str->histo, l_pressure_str->data_set, 
      g_pressure_str->min, g_pressure_str->max, l_indicies);
  compute_histo(l_density_str->histo, l_density_str->data_set, 
      g_density_str->min, g_density_str->max, l_indicies);
  compute_histo(l_temp_str->histo, l_temp_str->data_set, 
      g_temp_str->min, g_temp_str->max, l_indicies); 
  
  //Need to compute global histogram with local data in order to properly
  //allocate space for MPI_reduce
  //If space becomes an issue, we can limit populating global histo to only
  //process 0. We can then do an MPI_Reduce using local histos instead of global
  //if(rank==0){
    compute_histo(g_pressure_str->histo, l_pressure_str->data_set, 
        g_pressure_str->min, g_pressure_str->max, l_indicies);
    compute_histo(g_density_str->histo, l_density_str->data_set, 
        g_density_str->min, g_density_str->max, l_indicies);
    compute_histo(g_temp_str->histo, l_temp_str->data_set, 
        g_temp_str->min, g_temp_str->max, l_indicies);
  //}
  t4=MPI_Wtime();
  MPI_Allreduce(l_pressure_str->histo->count, g_pressure_str->histo->count, 
      g_pressure_str->histo->numcols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  t5=MPI_Wtime();
  MPI_Allreduce(l_temp_str->histo->count, g_temp_str->histo->count, 
      g_temp_str->histo->numcols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(l_density_str->histo->count, g_density_str->histo->count, 
      g_density_str->histo->numcols, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  t6=MPI_Wtime();

  if(rank==0){
    double time1, time2, time3, time4;
    time1=t2-t1;
    time2=(t3-t2)+time1;
    time3 = t5-t4;
    time4= (t6-t5)+time3;
    printf("Single Reduce All: %E\n",time1);
    printf("Nine Reduce All: %E\n",time2);
    printf("1 to 9 Scales: %f \n",time2/time1);
    printf("Single Reduce All Across Array: %E\n", time3);
    printf("3 Reduce All Across 3 Arrays: %E\n", time4);
    printf("1 to 3  Scales: %f \n",time4/time3);
  }
/*
  if(rank==0){
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
  }
*/

  //
  //********************************************************
  
  /////////////////////////////
  // Clean up heap variables
  /////////////////////////////

  //***************************************
  //Billy
  free_field_val(l_p);
  free_var_struct(l_pressure_str);
  free_var_struct(l_temp_str);
  free_var_struct(l_density_str);
  free_var_struct(g_pressure_str);
  free_var_struct(g_temp_str);
  free_var_struct(g_density_str);
  free_indicies(l_indicies);
      
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
