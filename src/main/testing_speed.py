
import time
import os
threads = 12
os.environ["MKL_NUM_THREADS"] = str(threads)
os.environ["NUMEXPR_NUM_THREADS"] = str(threads)
os.environ["OMP_NUM_THREADS"] = str(threads)
import numpy as np
# Generate random arrays A, B, C, and D

A = np.random.rand(100000, 50000)
B = np.random.rand(100000, 50000)
C = np.random.rand(100000)
D = np.random.rand(100000)

start_time = time.time()
result = np.sum(A* B, axis=1)/(C* D)
end_time = time.time()
print(result.shape)
print(end_time -start_time )


# Compute the result using einsum
start_time = time.time()
numerator = np.einsum('ij,ij->i', A, B)
denominator = np.einsum('i,i->i', C, D)
result_1 = numerator / denominator
end_time = time.time()

print(result_1.shape)
print(end_time -start_time )



from einsumt import einsumt as einsum
from einsumt import bench_einsumt
start_time = time.time()
numerator = bench_einsumt('ij,ij->i', A, B, pool =12)
denominator = bench_einsumt('i,i->i', C, D, pool= 2)
result_2 = numerator / denominator
end_time = time.time()

print(result_2.shape)
print(end_time -start_time)