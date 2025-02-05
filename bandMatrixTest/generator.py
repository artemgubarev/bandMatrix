import numpy as np

def generate_banded_matrix(n, b, low=1, high=10):
    random_values = np.random.randint(low, high, size=(n, 2 * b + 1))
    matrix = np.zeros((n, n))
    for i in range(n):
        start = max(0, i - b)
        end = min(n, i + b + 1)
        matrix[i, start:end] = random_values[i, :end - start]

    return matrix 

def generate_random_array(length, min_value=0, max_value=100):
    return [np.random.randint(min_value, max_value) for _ in range(length)]

n = 4000
b = 1756
banded_matrix = generate_banded_matrix(n, b)
c = generate_random_array(n, min_value=0.1, max_value=10)