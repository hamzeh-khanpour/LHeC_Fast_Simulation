import tensorflow as tf
tf.debugging.set_log_device_placement(True)

# Create a tensor and run on GPU
with tf.device('/GPU:0'):
    a = tf.constant([[1.0, 2.0], [3.0, 4.0]])
    b = tf.matmul(a, a)
print(b)
