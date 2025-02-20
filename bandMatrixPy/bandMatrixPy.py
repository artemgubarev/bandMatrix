#import subprocess

#num_processes = 4
#subprocess.run(["mpiexec", "-n", str(num_processes), "python", "run_mpi.py"])



#start_time = time.time()

#l,u = solver.lu_decomposition(A,b,n)
#x = solver.solve_lu(l,u,c,n)
##printer.print_1d(x, 'X',True)

#end_time = time.time()
#elapsed_time = end_time - start_time
#print(f"time: {elapsed_time:.6f} sec")



from fractions import Fraction
import numpy as np

def generate_random_array(length, min_value=0, max_value=100):
    return [np.random.randint(min_value, max_value) for _ in range(length)]

def write_matrix(filename, n, b, matrix, c):
    lines = [f"{n}\n", f"{b}\n"]
    lines.extend(' '.join(map(str, map(float, row))) + '\n' for row in matrix)
    lines.append(' '.join(map(str, map(float, c))) + '\n')
    with open(filename, 'w') as f:
        f.writelines(lines)

def generate_banded_matrix(n, b, low=1, high=10):
    random_values = np.random.randint(low, high, size=(n, 2 * b + 1))
    matrix = np.zeros((n, n))
    for i in range(n):
        start = max(0, i - b)
        end = min(n, i + b + 1)
        matrix[i, start:end] = random_values[i, :end - start]
    return matrix   

if __name__ == '__main__':
    n = 5000
    b = 2000
    banded_matrix = generate_banded_matrix(n, b)
    c = generate_random_array(n, min_value=0.1, max_value=10)
    write_matrix('matrix5000.txt', n,b,banded_matrix, c)