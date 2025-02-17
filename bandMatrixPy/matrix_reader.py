import numpy as np

def read_matrix(filename):
    with open(filename, 'r') as f:
        n = int(f.readline().strip())
        b = int(f.readline().strip())
        
        matrix = np.array([list(map(float, f.readline().split())) for _ in range(n)], dtype=np.float64).flatten()
        c = np.array(list(map(float, f.readline().split())), dtype=np.float64).flatten()
        
        if matrix.shape[0] != n * n:
            raise ValueError(f"Matrix must be of size {n}x{n}, but got mismatch in size.")
        if c.shape[0] != n:
            raise ValueError(f"Vector c must have {n} elements, but got {c.shape[0]}")
        
        return n, b, matrix.reshape(n, n), c
