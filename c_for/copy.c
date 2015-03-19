/*
 * In-situ analyses 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <pnetcdf.h>
#include "insitustats.h"

int dim1, dim2, iter;
static int timestep=5;

/*
 * copyfld function
 * - copies data for further analyses
 *
 *	Input parameters:
 *		rank - rank of the process
 *		sim_ts - simulation time step
 *		fname - field name
 *
 */

void pdata(int rank, int sim_ts,char *fname, int num_dims, int *dims, double *field, MPI_Comm comm){
  int dim1, dim2, status, tglob;
  static int ncid,avg_id,tgrid, tstep, min_id, max_id, g_avg_id, g_max_id, g_min_id;
  char *filename;
  char *varname;
  const char* extension = ".nc";
  char buffer[16];
  FILE *fp;
  int size;
  int dimids[2];
  MPI_Offset count[2];
  MPI_Offset index_v[2];
  MPI_Offset glob;
  
  sprintf(buffer, "%i", sim_ts);
  filename = malloc(strlen(fname)+strlen(extension)+4+1);
  strcpy(filename, fname);
  strcat(filename, buffer);
  strcat(filename, extension);

  field_val* g_temp = field_val_new_empty("GlobalTemperature");

  field_val* l_temp = field_val_new("LocalTemperature", field);
  
  compute_var_stats(l_temp, dims);

  MPI_Allreduce(get_mean_ptr(l_temp), get_mean_ptr(g_temp), 1,MPI_DOUBLE, 
      MPI_SUM, comm);
  MPI_Allreduce(get_min_ptr(l_temp), get_min_ptr(g_temp), 1,MPI_DOUBLE, 
      MPI_MIN, comm);
  MPI_Allreduce(get_max_ptr(l_temp), get_max_ptr(g_temp), 1,MPI_DOUBLE, 
      MPI_MAX, comm);
  MPI_Comm_size(comm, &size);
  *get_mean_ptr(g_temp)/=size;

  status = ncmpi_create(comm, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid);

  printf("%d:%d\n", rank, ncid);
  //status = ncmpi_def_dim (ncid, "TimeSteps",NC_UNLIMITED , &tstep);
  //status = ncmpi_def_dim (ncid, "TimeSteps",NC_UNLIMITED , &tstep);
  status = ncmpi_def_dim (ncid, "NumSubGrids", size, &tgrid); 
  status = ncmpi_def_dim (ncid, "Global Dimensions", 1, &tglob);

  dimids[0]=tstep;
  dimids[1]=tgrid;

  sprintf(buffer,"%i4", rank);
  varname = malloc(strlen(fname)+4+3+1);
  strcpy(varname, fname);
  strcat(varname, "avg");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tgrid,&avg_id);

  varname = realloc(varname, strlen(fname)+4+3+1);
  strcpy(varname,fname);
  strcat(varname, "min");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tgrid ,&min_id);

  varname = realloc(varname, strlen(fname)+4+3+1);
  strcpy(varname,fname);
  strcat(varname, "max");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tgrid ,&max_id);
  
  sprintf(buffer,"%i4", rank);
  varname = malloc(strlen(fname)+4+3+1);
  strcpy(varname, fname);
  strcat(varname, "Gavg");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tglob,&g_avg_id);

  varname = realloc(varname, strlen(fname)+4+3+1);
  strcpy(varname,fname);
  strcat(varname, "Gmin");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tglob ,&g_min_id);

  varname = realloc(varname, strlen(fname)+4+3+1);
  strcpy(varname,fname);
  strcat(varname, "Gmax");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tglob ,&g_max_id);
  
  free(varname);  

  status = ncmpi_enddef(ncid);
  
  index_v[0]=sim_ts;
  index_v[1]=rank;
  count[0]=1;
  count[1]=1;
  glob=0;

  //status = ncmpi_begin_indep_data(ncid); 
  //status = ncmpi_put_var1(ncid, avg_id, index_v, get_mean_ptr(l_temp) , count[1], MPI_DOUBLE);
  //status = ncmpi_end_indep_data(ncid);
  /*
  status = ncmpi_put_vara_all(ncid, avg_id, index_v, count, get_mean_ptr(l_temp), count[0], MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, min_id, index_v, count, get_min_ptr(l_temp), count[0], MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, max_id, index_v, count, get_max_ptr(l_temp), count[0], MPI_DOUBLE);
*/
  status = ncmpi_put_vara_all(ncid, avg_id, &index_v[1], &count[0], get_mean_ptr(l_temp), count[0], MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, min_id, &index_v[1], &count[0], get_min_ptr(l_temp), count[0], MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, max_id, &index_v[1], &count[0], get_max_ptr(l_temp), count[0], MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, g_avg_id, &glob, &count[0], get_mean_ptr(g_temp), count[0], MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, g_min_id, &glob, &count[0], get_min_ptr(g_temp), count[0], MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, g_max_id, &glob, &count[0], get_max_ptr(g_temp), count[0], MPI_DOUBLE);
  //status = ncmpi_sync_numrecs(ncid);
  
  //if(sim_ts ==  19) 
    
    status = ncmpi_close(ncid);


  printf("Rank:%d L Mean:%f L Min:%f L Max:%f G Mean:%f G Min:%f G Max:%f\n", rank, *get_mean_ptr(l_temp), 
      *get_min_ptr(l_temp), *get_max_ptr(l_temp),  *get_mean_ptr(g_temp), *get_min_ptr(g_temp), *get_max_ptr(g_temp));

  free_field_val(g_temp);
  free_field_val(l_temp);

  //fp = fopen(filename, "w");
  
  if (rank < 20) { 
    for(dim2=0; dim2<dims[2] ; dim2++) { 
      for(dim1=0; dim1<dims[1] ; dim1++) {
      
        int iter = dim1 + dim2 * (dims[1]-dims[0]+1);
	//fprintf(fp, "%lf\n", field[iter]);
        //fprintf(fp, "%lf\n",field[iter]);
      }
    }
  }

  //free(filename);
  //fclose(fp);

}

//void copyfld(int rank, int timestep, char *fname, int beg_dim1, int end_dim1, int beg_dim2, int end_dim2, double *field) {
void cpdata(int rank, int sim_ts, char *fname, int num_dims, int *dims, double *field, MPI_Fint* com_group) {

  int dim1, dim2;
  int ranktest;
  MPI_Comm c_com_group;
  int size;
  c_com_group = MPI_Comm_f2c(*com_group);
  /*if (timestep == sim_ts) {
    printf ("Timestep %d sim_ts %d\n", timestep, sim_ts);
    return;
  }
  else
    timestep = sim_ts;
  */
  if (rank < 20){ 
    //printf("%d: %d copy %s %d %d %d %d\n", rank, timestep, fname, dims[0], dims[1], dims[2]);
  }
  //printf("%d\n",dims[2]);

  //printf("%d\n",dims[1]);

  if (rank < 20) { 
    for(dim2=0; dim2<dims[1] ; dim2++) { 
      for(dim1=0; dim1<dims[0] ; dim1++) {
			//int iter = dim1 + dim2 * (dims[1]-dims[0]+1);
                        int iter = dim2*dims[1] +dim1;
                        //printf("%d\n", iter);
			//printf("%d: %d val %s[%d] %lf\n", rank, sim_ts, fname, iter, field[iter]);
      }
    }
  }
  if(rank ==0 )
    ranktest=rank;

  MPI_Bcast(&ranktest, 1,MPI_INTEGER, 0, c_com_group);

  //printf("%i:%i\n", ranktest, rank);

  pdata(rank, sim_ts, fname, num_dims, dims, field, c_com_group);
} 



//called from phys_run1
void copyfield(int *sim_ts, double *ps, double *t, double *u, double *v) {
	
	printf("copyfield: %d %lf %lf %lf %lf\n", *sim_ts, ps[0], t[0], u[0], v[0]);

	return;

}

