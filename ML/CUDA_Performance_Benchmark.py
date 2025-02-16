import tensorflow as tf
import time

# Generate random data
x = tf.random.normal([10000, 10000])

# Time CPU execution
start_cpu = time.time()
_ = tf.matmul(x, x)
end_cpu = time.time()
print(f"CPU computation time: {end_cpu - start_cpu:.4f} seconds")

# Time GPU execution
with tf.device('/GPU:0'):
    start_gpu = time.time()
    _ = tf.matmul(x, x)
    end_gpu = time.time()
print(f"GPU computation time: {end_gpu - start_gpu:.4f} seconds")
