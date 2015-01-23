# 	---------------------------
# 	Compiler Flags
# 	---------------------------
MPICC =	mpicc
MPICXX = mpicxx

#Compiler Flags
MPI_LDFLAGS = 
MPI_CFLAGS += -g -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE 

CXXFLAGS = $(MPI_CFLAGS)
LDFLAGS = $(MPI_LDFLAGS)

all: 3D_Grid 

%.o: %.c 
	$(MPICC) $(MPI_CFLAGS) $< -c -o $@

3D_Grid: parallel_3D.o 
	$(MPICC) -o 3D_Grid parallel_3D.o -lm	

clean:
	rm -f *.a *.o a.out core* 3D_Grid
