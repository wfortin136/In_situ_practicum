# ---------------------------
# 	Compiler Flags
# 	---------------------------
MPICC =	mpicc
MPICXX = mpicxx
MPIF = mpifort

#Compiler Flags
MPI_LDFLAGS = 
MPI_CFLAGS += -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE 

CXXFLAGS = $(MPI_CFLAGS)
LDFLAGS = $(MPI_LDFLAGS)

all: run 
%.o: %.f90
	$(MPIF) $< -c -o $@

%.o: %.c 
	$(MPICC) $(MPI_CFLAGS) $< -c -I/Users/billyfortin/PnetCDF/include/ -o $@

run: insitu_for.o copy.o insitustats.o 
	$(MPIF) -o run insitu_for.o copy.o insitustats.o -L/Users/billyfortin/PnetCDF/lib -lpnetcdf	

clean:
	rm -f *.a *.o a.out core* run
