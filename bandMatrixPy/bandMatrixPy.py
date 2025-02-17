import subprocess

num_processes = 4
subprocess.run(["mpiexec", "-n", str(num_processes), "python", "run_mpi.py"])



#start_time = time.time()

#l,u = solver.lu_decomposition(A,b,n)
#x = solver.solve_lu(l,u,c,n)
##printer.print_1d(x, 'X',True)

#end_time = time.time()
#elapsed_time = end_time - start_time
#print(f"time: {elapsed_time:.6f} sec")