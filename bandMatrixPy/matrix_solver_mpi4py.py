import numpy as np
from mpi4py import MPI

def lu_decomposition_banded(m, b, n, comm=MPI.COMM_WORLD):
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    rows_per_proc = n // size
    if n % size != 0:
        raise ValueError("n must be divisible by the number of processes.")
    
    start_row = rank * rows_per_proc
    end_row = (rank + 1) * rows_per_proc
    
    L_local = np.zeros_like(m)
    U_local = np.copy(m)
    
    for k in range(n - 1):
        owner_proc = k // rows_per_proc
        
        if rank == owner_proc:
            local_k = k - start_row
            pivot_row = U_local[local_k, :].copy()
        else:
            pivot_row = np.zeros(n, dtype=m.dtype)
        
        comm.Bcast(pivot_row, root=owner_proc)
        
        i_min = max(k + 1, start_row)
        i_max = min(k + b + 1, end_row)
        
        for i in range(i_min, i_max):
            local_i = i - start_row
            
            if abs(pivot_row[k]) < 1e-15:
                raise ZeroDivisionError("Zero pivot encountered.")
            
            L_local[local_i, k] = U_local[local_i, k] / pivot_row[k]
            j_max = min(k + b + 1, n)
            factor = L_local[local_i, k]
            U_local[local_i, k:j_max] -= factor * pivot_row[k:j_max]
        
        comm.Barrier()
    
    for i_loc in range(L_local.shape[0]):
        global_i = start_row + i_loc
        if global_i < n:
            L_local[i_loc, global_i] = 1.0
    
    return L_local, U_local, start_row

def forward_substitution_banded(L_local, b_vec_local, n, comm=MPI.COMM_WORLD):
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    rows_per_proc = n // size
    start_row = rank * rows_per_proc
    end_row = (rank + 1) * rows_per_proc
    
    y_local = np.zeros_like(b_vec_local)
    
    for i in range(n):
        owner_proc = i // rows_per_proc
        
        if rank == owner_proc:
            local_i = i - start_row
            sum_part = 0.0
            
            for j in range(i):
                owner_proc_j = j // rows_per_proc
                if owner_proc_j == rank:
                    local_j = j - start_row
                    sum_part += L_local[local_i, j] * y_local[local_j]
                else:
                    y_j = np.array([0.0], dtype=b_vec_local.dtype)
                    comm.Recv(y_j, source=owner_proc_j, tag=j)
                    sum_part += L_local[local_i, j] * y_j[0]
            
            y_local[local_i] = b_vec_local[local_i] - sum_part
            y_to_send = np.array([y_local[local_i]], dtype=b_vec_local.dtype)
            
            for dest in range(size):
                if dest != rank:
                    comm.Send(y_to_send, dest=dest, tag=i)
        else:
            y_i = np.array([0.0], dtype=b_vec_local.dtype)
            comm.Recv(y_i, source=owner_proc, tag=i)
            for dest in range(size):
                if dest != owner_proc:
                    comm.Send(y_i, dest=dest, tag=i)
    
    return y_local

def backward_substitution_banded(U_local, y_vec_local, n, comm=MPI.COMM_WORLD):
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    rows_per_proc = n // size
    start_row = rank * rows_per_proc
    end_row = (rank + 1) * rows_per_proc
    
    x_local = np.zeros_like(y_vec_local)
    
    for i in reversed(range(n)):
        owner_proc = i // rows_per_proc
        
        if rank == owner_proc:
            local_i = i - start_row
            sum_part = 0.0
            
            for j in range(i + 1, n):
                owner_proc_j = j // rows_per_proc
                if owner_proc_j == rank:
                    local_j = j - start_row
                    sum_part += U_local[local_i, j] * x_local[local_j]
                else:
                    x_j = np.array([0.0], dtype=y_vec_local.dtype)
                    comm.Recv(x_j, source=owner_proc_j, tag=n + j)
                    sum_part += U_local[local_i, j] * x_j[0]
            
            if abs(U_local[local_i, i]) < 1e-15:
                raise ZeroDivisionError("Zero pivot encountered.")
            
            x_local[local_i] = (y_vec_local[local_i] - sum_part) / U_local[local_i, i]
            x_to_send = np.array([x_local[local_i]], dtype=y_vec_local.dtype)
            
            for dest in range(size):
                if dest != rank:
                    comm.Send(x_to_send, dest=dest, tag=n + i)
        else:
            x_i = np.array([0.0], dtype=y_vec_local.dtype)
            comm.Recv(x_i, source=owner_proc, tag=n + i)
            for dest in range(size):
                if dest != owner_proc:
                    comm.Send(x_i, dest=dest, tag=n + i)
    
    return x_local

def solve_system_banded(m_global, b_global, n, bandwidth):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    
    if m_global.shape != (n, n):
        raise ValueError("m_global must be n x n.")
    if b_global.shape != (n,):
        raise ValueError("b_global must be n-dimensional.")
    
    rows_per_proc = n // size
    start_row = rank * rows_per_proc
    end_row = (rank + 1) * rows_per_proc
    
    m_local = m_global[start_row:end_row, :].copy()
    b_local = b_global[start_row:end_row].copy()
    
    L_local, U_local, _ = lu_decomposition_banded(m_local, bandwidth, n, comm)
    y_local = forward_substitution_banded(L_local, b_local, n, comm)
    x_local = backward_substitution_banded(U_local, y_local, n, comm)
    
    x_gathered = None
    if rank == 0:
        x_gathered = np.zeros(n, dtype=x_local.dtype)
    
    comm.Gather(x_local, x_gathered, root=0)
    
    return x_gathered