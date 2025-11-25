#ifndef IMAGE_COMPRESSION_H
#define IMAGE_COMPRESSION_H

/*
 * JPEG Compression Library
 * 
 * This library implements basic JPEG compression and decompression.
 * It supports encoding and decoding of RGB images using:
 * - DCT (Discrete Cosine Transform)
 * - Standard JPEG quantization matrices
 * - Huffman encoding with standard JPEG tables
 * - YCbCr color space conversion
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>

/* Utility macros */
#define min(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a < _b ? _a : _b; })

#define max(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })

/* ============================================================================
 * CONSTANTS - Quantization matrices and Huffman tables
 * ============================================================================ */

extern const float LUMINANCE_QUANT[64];
extern const float CHROMINANCE_QUANT[64];
extern const uint8_t ZIGZAG_IDX[64];
extern const uint16_t AC_LUMINANCE_CODES[162];
extern const uint8_t AC_LUMINANCE_LENGTHS[162];
extern const uint16_t AC_CHROMINANCE_CODES[162];
extern const uint8_t AC_CHROMINANCE_LENGTHS[162];
extern const uint8_t AC_VALUES[162];
extern const uint16_t DC_LUMINANCE_CODES[12];
extern const uint8_t DC_LUMINANCE_LENGTHS[12];
extern const uint16_t DC_CHROMINANCE_CODES[12];
extern const uint8_t DC_CHROMINANCE_LENGTHS[12];
extern const uint8_t DC_VALUES[12];

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================ */

/* Internal structure for entropy encoding */
typedef struct DATA_PACKET {
    uint8_t rrrrssss;
    uint16_t rs_code;
    uint16_t VAL;
    int VAL_bits, rs_code_bits;
} DATA_PACKET;

/* Output bit stream for encoding */
typedef struct OUTSTREAM {
    const char* filename;
    FILE *file;
    int written_bits, written_bytes, buffer_bytes;
    char *buffer;
} OUTSTREAM;

/* Input bit stream for decoding */
typedef struct INSTREAM {
    const char* filename;
    FILE *file;
    int read_bits, read_bytes, buffer_bytes;
    char *buffer;
    bool eof;
} INSTREAM;

/* RGB image container */
typedef struct RGB_IMAGE {
    uint8_t *r, *g, *b;
    uint16_t HEIGHT, WIDTH;
} RGB_IMAGE;

/* ============================================================================
 * PUBLIC API - Main encoding/decoding functions
 * ============================================================================ */

/**
 * Encode an RGB image to a compressed file.
 * 
 * @param filename Output file path
 * @param r Red channel data (modified during encoding)
 * @param g Green channel data (modified during encoding)
 * @param b Blue channel data (modified during encoding)
 * @param width Image width
 * @param height Image height
 * @return 0 on success, -1 on error
 */
int encode_image(char *filename, uint8_t *r, uint8_t *g, uint8_t *b, uint16_t width, uint16_t height);

/**
 * Decode a compressed file to RGB image.
 * 
 * @param filename Input file path
 * @param r Pointer to receive red channel (allocated by function)
 * @param g Pointer to receive green channel (allocated by function)
 * @param b Pointer to receive blue channel (allocated by function)
 * @param width Pointer to receive image width
 * @param height Pointer to receive image height
 * @return 0 on success, -1 on error
 */
int decode_image(char *filename, uint8_t **r, uint8_t **g, uint8_t **b, uint16_t *width, uint16_t *height);

/* RGB image allocation/deallocation */
RGB_IMAGE *new_RGB_IMAGE(uint16_t width, uint16_t height);
void delete_RGB_IMAGE(RGB_IMAGE *img);

/* ============================================================================
 * INTERNAL API - Used by the compression pipeline
 * ============================================================================ */

/* Block extraction/placement */
void get_block(uint8_t *slice, uint8_t *UINT8_BLOCK, uint16_t IMG_WIDTH, uint16_t IMG_HEIGHT, int I0, int J0);
void put_block(uint8_t *slice, uint8_t *UINT8_BLOCK, uint16_t IMG_WIDTH, uint16_t IMG_HEIGHT, int I0, int J0);

/* Color space conversion */
void image_rgb_to_yCbCr(uint8_t *r_to_y, uint8_t *g_to_Cb, uint8_t *b_to_Cr);
void image_yCbCr_to_rgb(uint8_t *y_to_r, uint8_t *Cb_to_g, uint8_t *Cr_to_b);
void slice_rgb_to_yCbCr(uint8_t *r_to_y, uint8_t *g_to_Cb, uint8_t *b_to_Cr, size_t N);
void slice_yCbCr_to_rgb(uint8_t *y_to_r, uint8_t *Cb_to_g, uint8_t *Cr_to_b, size_t N);
void downsample_420(uint8_t *Cb, uint8_t *Cr, size_t N);

/* DCT operations */
void block_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK);
void block_inv_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK);
void general_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_WIDTH, size_t BLOCK_HEIGHT);

/* Quantization */
void block_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK);
void block_inv_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK);

/* Zigzag serialization */
void block_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX);
void block_inv_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX);

/* Entropy coding */
void DATA_PACKET_pack(DATA_PACKET *data, int16_t VAL, uint8_t zeros);
bool DATA_PACKET_encode(DATA_PACKET *data, const uint16_t *CODES, const uint8_t *VALUES, const uint8_t *LENGTHS, size_t N_CODES);
int min_bits_abs(int16_t n);

int block_encode(OUTSTREAM* out, int16_t *INT16_SEQUENCE, int16_t *PREV_DC,
    const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS,
    const uint16_t *AC_CODES, const uint8_t *AC_VALUES, const uint8_t *AC_LENGTHS);

int block_decode(INSTREAM* in, int16_t *INT16_SEQUENCE, int16_t* PREV_DC,
    const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS,
    const uint16_t *AC_CODES, const uint8_t *AC_VALUES, const uint8_t *AC_LENGTHS);

bool search_codes(INSTREAM *in, uint8_t *rrrrssss, const uint16_t *CODES, const uint8_t *VALUES, const uint8_t *LENGTHS, size_t CODES_NUMBER);
void decode_data(INSTREAM *in, uint8_t rrrrssss, int *ssss, int *rrrr, uint16_t *val);
bool write_data(int16_t *INT16_SEQUENCE, bool is_dc, int idx, int ssss, int rrrr, uint16_t val);

/* Slice-level encoding/decoding */
int encode_slice(OUTSTREAM *out, uint8_t *slice, uint16_t WIDTH, uint16_t HEIGHT, bool is_luminance);
int decode_slice(INSTREAM *in, uint8_t *slice, uint16_t WIDTH, uint16_t HEIGHT, bool is_luminance);

/* Bit stream I/O */
OUTSTREAM *new_OUTSTREAM(const char* filename, int buffer_bytes);
int delete_OUTSTREAM(OUTSTREAM *out);
int OUTSTREAM_push(OUTSTREAM *out, uint16_t data, int bits);
void OUTSTREAM_reset_bytes(OUTSTREAM *out);

INSTREAM *new_INSTREAM(const char* filename, int buffer_bytes);
int delete_INSTREAM(INSTREAM *in);
int INSTREAM_pull(INSTREAM *in, uint16_t *data, int bits);
int INSTREAM_pull_1bit(INSTREAM *in);
void INSTREAM_reset_bytes(INSTREAM *in);

#endif /* IMAGE_COMPRESSION_H */