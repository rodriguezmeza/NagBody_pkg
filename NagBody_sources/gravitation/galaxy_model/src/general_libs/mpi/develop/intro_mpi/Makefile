# Makefile to build mpi class examples
# --- macros
CC= mpicc
FC=mpif90

# --- targets
all:    hello_world_c hello_world_f numerical_integ numerical_integ2 matmult array_search1 \
array_search2 numerical_integ_err

hello_world_c:   hello_world_c.c
	$(CC)  -o hello_world_c  hello_world_c.c
	
hello_world_f:   hello_world_f.f
	$(FC)  -o hello_world_f  hello_world_f.f
	
numerical_integ: numerical_integ.c
	$(CC) -o numerical_integ  numerical_integ.c
	
numerical_integ2: numerical_integ2.c
	$(CC) -o numerical_integ2  numerical_integ2.c

matmult: matmult.f
	$(FC)  -o matmult  matmult.f
	
array_search1: array_search1.c
	$(CC)  -o array_search1  array_search1.c
	
array_search2: array_search2.c
	$(CC)  -o array_search2  array_search2.c
 
numerical_integ_err: numerical_integ_err.c
	$(CC)  -o numerical_integ_err  numerical_integ_err.c
	
# --- remove binary and object files
clean:
	rm -f *.o hello_world_c hello_world_f numerical_integ numerical_integ2 matmult array_search1 \
array_search2 numerical_integ_err
