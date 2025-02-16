import tensorflow as tf
import time

# Check GPU availability
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))

# Test CUDA performance
start_time = time.time()
with tf.device('/GPU:0'):
    a = tf.random.normal([10000, 10000])
    b = tf.random.normal([10000, 10000])
    c = tf.matmul(a, b)

end_time = time.time()
print("GPU computation time:", end_time - start_time, "seconds")
