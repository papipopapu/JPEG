from shitpeg import *
import numpy as np

# create 3 random numpy arrays: r,g and b, of shape (8,8)
r = np.random.randint(0, 255, (8,8), dtype = np.uint8)
g = np.random.randint(0, 255, (8,8), dtype = np.uint8)
b = np.random.randint(0, 255, (8,8), dtype = np.uint8)



rgb_in = [r, g, b]

ENCODE_RGB("pyc.bin", rgb_in)
print(r)
print(g)
print(b)
#rgb_out = DECODE_RGB("test.bin")

