# JPEG Compression Library

A C implementation of the JPEG compression algorithm, supporting encoding and decoding of RGB images using the standard JPEG pipeline.

## Features

- **DCT (Discrete Cosine Transform)**: 8x8 block-based DCT and inverse DCT
- **Quantization**: Standard JPEG luminance and chrominance quantization matrices
- **Huffman Coding**: Standard JPEG Huffman tables for AC and DC coefficients
- **Color Space Conversion**: RGB to YCbCr and back
- **4:2:0 Chroma Subsampling**: Reduces chrominance data for better compression
- **Zigzag Serialization**: Standard JPEG zigzag pattern for coefficient ordering

## Project Structure

```
JPEG/
├── CMakeLists.txt          # CMake build configuration
├── README.md               # This file
├── main.c                  # Example/test program
├── core/                   # Core library sources
│   ├── image_compression.h # Main header file
│   ├── block_process.c     # Block-level operations (DCT, quantization, encoding)
│   ├── image_process.c     # Image-level operations (color conversion, slice processing)
│   ├── streams.c           # Bit stream I/O handling
│   └── constants.c         # Huffman tables and quantization matrices
```

## Building

### Prerequisites

- CMake 3.10 or later
- C99 compatible compiler (GCC, Clang, etc.)

### Build Instructions

```bash
# Create build directory
mkdir build && cd build

# Configure with CMake
cmake ..

# Build
make

# Run the test program
./jpeg_test
```

### Build Options

```bash
# Release build with optimizations
cmake -DCMAKE_BUILD_TYPE=Release ..

# Debug build
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

## Usage

### API Overview

The library provides the following main functions:

#### Image Encoding
```c
int encode_image(char *filename, uint8_t *r, uint8_t *g, uint8_t *b, 
                 uint16_t width, uint16_t height);
```
Encodes an RGB image to a compressed binary file.

#### Image Decoding
```c
int decode_image(char *filename, uint8_t **r, uint8_t **g, uint8_t **b,
                 uint16_t *width, uint16_t *height);
```
Decodes a compressed binary file back to RGB image data.

### Example

```c
#include "core/image_compression.h"

int main() {
    // Create sample image data (8x8 pixels)
    uint8_t *r = malloc(64);
    uint8_t *g = malloc(64);
    uint8_t *b = malloc(64);
    
    // ... fill with image data ...
    
    // Encode to file
    encode_image("output.bin", r, g, b, 8, 8);
    
    // Decode from file
    uint8_t *r_out, *g_out, *b_out;
    uint16_t width, height;
    decode_image("output.bin", &r_out, &g_out, &b_out, &width, &height);
    
    // Clean up
    free(r); free(g); free(b);
    free(r_out); free(g_out); free(b_out);
    
    return 0;
}
```

## JPEG Pipeline

The implementation follows the standard JPEG compression pipeline:

### Encoding
1. **Color Conversion**: RGB → YCbCr
2. **Chroma Subsampling**: 4:2:0 subsampling on Cb and Cr channels
3. **Block Division**: Split image into 8×8 blocks
4. **DCT**: Apply 2D Discrete Cosine Transform to each block
5. **Quantization**: Divide by quantization matrix and round
6. **Zigzag Scan**: Reorder coefficients in zigzag pattern
7. **Entropy Coding**: Huffman encode DC (differential) and AC coefficients

### Decoding
The decoding process reverses these steps.

## File Format

The compressed file format is a custom binary format (not standard JFIF/JPEG):
- Bytes 0-1: Image width (16-bit)
- Bytes 2-3: Image height (16-bit)
- Remaining: Huffman-encoded Y, Cb, Cr data

## Limitations

- Custom binary format (not compatible with standard JPEG viewers)
- No support for progressive encoding
- Fixed quantization quality
- Image dimensions should ideally be multiples of 8

## License

This project is provided as-is for educational purposes.
