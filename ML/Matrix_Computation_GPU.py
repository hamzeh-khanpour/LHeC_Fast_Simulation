import tensorflow as tf
import time

# Create large matrices
A = tf.random.normal([1000, 1000])
B = tf.random.normal([1000, 1000])

# Measure CPU time
start = time.time()
C = tf.matmul(A, B)
cpu_time = time.time() - start

# Measure GPU time
with tf.device('/GPU:0'):
    start = time.time()
    C_gpu = tf.matmul(A, B)
    gpu_time = time.time() - start

print(f"CPU computation time: {cpu_time:.4f} seconds")
print(f"GPU computation time: {gpu_time:.4f} seconds")
