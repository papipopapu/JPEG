from ctypes import CDLL, POINTER, c_int, c_char_p
from numpy import ndarray, uint8, uint16, ctypeslib, asarray
import os

SHITPEG = CDLL(f"{os.path.dirname(__file__)}/libshitpeg.so")
PTR8 = ctypeslib.ndpointer(dtype = uint8, ndim = 2, flags = "C")
PTR16 = ctypeslib.ndpointer(dtype = uint16, ndim = 1, flags = "C")


SHITPEG.encode_image.argtypes = [c_char_p, PTR8, PTR8, PTR8, PTR16]
SHITPEG.encode_image.restype = c_int

SHITPEG.decode_image.argtypes = [c_char_p, PTR8, PTR8, PTR8, PTR16]
SHITPEG.decode_image.restype = c_int

def ENCODE_RGB(filename, rgb):   
    SHITPEG.encode_image(filename.encode("utf-8"), rgb[0], rgb[1], rgb[2], asarray(rgb[0].shape, dtype = uint16))

def DECODE_RGB(filename):
    shape = ndarray(shape = (2,), dtype = uint16)
    r = ndarray(shape = (shape[0],shape[1]), dtype = uint8)
    g = ndarray(shape = (shape[0],shape[1]), dtype = uint8)
    b = ndarray(shape = (shape[0],shape[1]), dtype = uint8)

    SHITPEG.decode_image(filename.encode("utf-8"), r, g, b, shape)
    return [r, g, b]



