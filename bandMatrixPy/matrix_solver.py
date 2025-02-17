import numpy as np

def lu_decomposition(m, b, n):
    u = m.copy()
    l = np.eye(n, dtype=np.float64)
    
    for k in range(n - 1):
        for i in range(k + 1, min(k + b + 1, n)):
            l[i, k] = u[i, k] / u[k, k]
            u[i, k:min(k + b + 1, n)] -= l[i, k] * u[k, k:min(k + b + 1, n)]
    
    return l, u

def solve_lu(l, u, c, n):
    y = np.zeros(n, dtype=np.float64)
    for i in range(n):
        y[i] = c[i] - np.dot(l[i, :i], y[:i])
    
    x = np.zeros(n, dtype=np.float64)
    for i in reversed(range(n)):
        x[i] = (y[i] - np.dot(u[i, i + 1:], x[i + 1:])) / u[i, i]
    
    return x
