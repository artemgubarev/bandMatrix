import numpy as np
import sys
import os

def generate_banded_matrix(n, b, low=1, high=10):
    random_values = np.random.randint(low, high, size=(n, 2 * b + 1))
    matrix = np.zeros((n, n))
    for i in range(n):
        start = max(0, i - b)
        end = min(n, i + b + 1)
        matrix[i, start:end] = random_values[i, :end - start]
    return matrix 

def generate_random_array(length, min_value=0, max_value=100):
    return [np.random.uniform(min_value, max_value) for _ in range(length)]

def write_matrix(filename, n, b, matrix, c):
    lines = [f"{n}\n", f"{b}\n"]
    lines.extend(' '.join(map(str, row)) + '\n' for row in matrix)
    lines.append(' '.join(map(str, c)) + '\n')
    with open(filename, 'w') as f:
        f.writelines(lines)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 generate_matrix.py <n> <b>")
        sys.exit(1)

    n = int(sys.argv[1])
    b = int(sys.argv[2])

    filename = f"matrix{n}x{b}.txt"
    filepath = os.path.join("testData", filename)
    matrix = generate_banded_matrix(n, b)
    c = generate_random_array(n, min_value=0.1, max_value=10)
    write_matrix(filepath, n, b, matrix, c)