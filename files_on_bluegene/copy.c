/*
 * In-situ analyses 
 */

#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "insitustats.h"
#include <pnetcdf.h>

int dim1, dim2, iter;
static int timestep=0;

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
void pdata(int rank, int sim_ts,char *fname, int num_dims, int *dims, double *field, MPI_Comm com){
  int dim1, dim2, size, status,  tglob;
  static int ncid, avg_id, tstat, tstep, min_id, max_id, g_avg_id, g_max_id, g_min_id;
  char *filename;
  char *varname;
  const char* path = "/projects/ExaHDF5/fortin/archive/test123/atm/hist/";
  const char* extension = ".nc";
  char buffer[16];
  int* l_indicies;
  MPI_Offset start, count=1;
  MPI_Offset glob;
 /*
  sprintf(buffer, "%i", sim_ts);
  filename = malloc(strlen(fname)+strlen(extension)+4+1);
  strcpy(filename, fname);
  strcat(filename, buffer);
  strcat(filename, extension);
*/
  field_val* g_temp = field_val_new_empty("Global Temperature");
  field_val* l_temp = field_val_new("Local Temperature", field);	
  l_indicies = set_indicies(dims[1],dims[3], 1);	
  compute_var_stats(l_temp, l_indicies);
  MPI_Allreduce(get_mean_ptr(l_temp), get_mean_ptr(g_temp), 1,MPI_DOUBLE, 
	MPI_SUM, com);
  MPI_Allreduce(get_min_ptr(l_temp), get_min_ptr(g_temp), 1,MPI_DOUBLE, 
	MPI_MIN, com);
  MPI_Allreduce(get_max_ptr(l_temp), get_max_ptr(g_temp), 1,MPI_DOUBLE, 
	MPI_MAX, com);
  MPI_Comm_size(com, &size);
  *get_mean_ptr(g_temp)/=size;
 /*
  status = ncmpi_create(com, filename, NC_CLOBBER, MPI_INFO_NULL, &ncid); 

  status = ncmpi_def_dim (ncid, "NumSubGrids", size, &tstat);
  status = ncmpi_def_dim (ncid, "Global_Dimensions", 1, &tglob);

  varname = malloc(strlen(fname)+3+1);
  strcpy(varname, fname);
  strcat(varname, "avg");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tstat ,&avg_id);

  varname = (char*)realloc(varname, strlen(fname)+3+1);
  strcpy(varname,fname);
  strcat(varname, "min");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tstat ,&min_id);

  varname = (char*)realloc(varname, strlen(fname)+3+1);
  strcpy(varname,fname);
  strcat(varname, "max");
  status = ncmpi_def_var(ncid, varname, NC_DOUBLE, 1, &tstat ,&max_id);

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
  
  start=rank;
  count = 1;
  glob=0;

  status = ncmpi_put_vara_all(ncid, avg_id, &start, &count, get_mean_ptr(l_temp), count, MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, min_id, &start, &count, get_min_ptr(l_temp), count, MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, max_id, &start, &count, get_max_ptr(l_temp), count, MPI_DOUBLE);

  status = ncmpi_put_vara_all(ncid, g_avg_id, &glob, &count, get_mean_ptr(g_temp), count, MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, g_min_id, &glob, &count, get_min_ptr(g_temp), count, MPI_DOUBLE);
  status = ncmpi_put_vara_all(ncid, g_max_id, &glob, &count, get_max_ptr(g_temp), count, MPI_DOUBLE);
  
  status = ncmpi_close(ncid);
*/
  /*fp = fopen(filename, "w");
 
  fprintf(fp,"Rank:%d L Mean:%f L Min:%f L Max:%f\n", rank, *get_mean_ptr(l_temp), *get_min_ptr(l_temp), *get_max_ptr(l_temp));
  fprintf(fp,"Rank:%d G Mean:%f G Min:%f G Max:%f\n", rank, *get_mean_ptr(g_temp), *get_min_ptr(g_temp), *get_max_ptr(g_temp));
   
  if (rank < 20) { 
    for(dim2=0; dim2<dims[3] ; dim2++) { 
      for(dim1=0; dim1<dims[1] ; dim1++) {
        int iter = dim1 + dim2 * (dims[1]-dims[0]+1);
	//fprintf(fp, "%d: %d val %s[%d] %lf\n", rank, sim_ts, fname, iter, field[iter]);
      }
    }
  }
*/
  free_field_val(g_temp);
  free_field_val(l_temp);

  //free(filename);
  //fclose(fp);
}

//void copyfld(int rank, int timestep, char *fname, int beg_dim1, int end_dim1, int beg_dim2, int end_dim2, double *field) {
void cpdata(int rank, int sim_ts, char *fname, int num_dims, int *dims, double *field,  MPI_Fint* com_group) {

	int dim1, dim2;
	int * l_indicies;
  MPI_Comm c_com_group;
  c_com_group = MPI_Comm_f2c(*com_group);

  if (timestep == sim_ts) {
		//printf ("Timestep %d sim_ts %d\n", timestep, sim_ts);
		return;
  }
  else
		timestep = sim_ts;
/*
  if (rank < 20) 
		printf("%d: %d copy %s %d %d %d %d\n", rank, timestep, fname, dims[0], dims[1], dims[2], dims[3]);

  if (rank < 20) { 
	  for(dim2=0; dim2<dims[3] ; dim2++) { 
	   for(dim1=0; dim1<dims[1] ; dim1++) {
			int iter = dim1 + dim2 * (dims[1]-dims[0]+1);
			printf("%d: %d val %s[%d] %lf\n", rank, sim_ts, fname, iter, field[iter]);
		 }
	  }
	}
*/	
	pdata(rank, sim_ts, fname, num_dims, dims, field, c_com_group);
} 



//called from phys_run1
void copyfield(int *sim_ts, double *ps, double *t, double *u, double *v) {
	
	printf("copyfield: %d %lf %lf %lf %lf\n", *sim_ts, ps[0], t[0], u[0], v[0]);

	return;

}

