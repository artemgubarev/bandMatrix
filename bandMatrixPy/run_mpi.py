import matrix_reader as reader
import matrix_solver_mpi4py as mpi_solver
import time
from mpi4py import MPI

filename = 'matrix100.txt'
n,b,A,c = reader.read_matrix(filename)

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
    
if rank == 0:
    A_flat = A.flatten()
else:
    A_flat = None

if rank == 0:
    c_flat = c
else:
    c_flat = None

A_flat = comm.bcast(A_flat, root=0)
c_flat = comm.bcast(c_flat, root=0)

if A_flat is not None:
    A_global = A_flat.reshape(n, n)
else:
    A_global = None

if c_flat is not None:
    c_global = c_flat
else:
    c_global = None
    

start_time = time.time()

x_solution = mpi_solver.solve_system_banded(A_global, c_global, n, b)
#print(x_solution)

end_time = time.time()
elapsed_time = end_time - start_time
print(f"time: {elapsed_time:.6f} sec")
