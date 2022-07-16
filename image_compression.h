#ifndef IMAGE_COMPRESSION_H
#define IMAGE_COMPRESSION_H

// Includes
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>

#define min(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a < _b ? _a : _b; })

#define max(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })
// Constants
// BLOCK DIMENSIONS = 8X8
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
void print_matrix(int16_t *seq);
// data_node
typedef struct DATA_PACKET {
    /* Temporary struct to pack both AC and DC components */
    uint8_t rrrrssss; // packed here
    uint16_t rs_code;
    uint16_t VAL; // prob less than 16 bits, min bits will be packed and then recasted to an int16 to be interpreted
    int VAL_bits, rs_code_bits;

} DATA_PACKET;
typedef struct OUTSTREAM {
    const char* filename;
    FILE *file;
    int written_bits, written_bytes, buffer_bytes;
    // written means alredy written ~ dtr
    char *buffer;
} OUTSTREAM;
typedef struct INSTREAM {
    const char* filename;
    FILE *file;
    int read_bits, read_bytes, buffer_bytes;
    // written means alredy written ~ dtr
    char *buffer;
    bool eof;
} INSTREAM;

typedef struct RGB_IMAGE {
    uint8_t *r, *g, *b;
    uint16_t HEIGHT, WIDTH;
} RGB_IMAGE;
// block_process
void get_block(uint8_t *slice, uint8_t *UINT8_BLOCK, uint16_t IMG_WIDTH, uint16_t IMG_HEIGHT, int I0, int J0);
void put_block(uint8_t *slice, uint8_t *UINT8_BLOCK, uint16_t IMG_WIDTH, uint16_t IMG_HEIGHT, int I0, int J0);

void print_block(uint8_t *UINT8_BLOCK);
void image_rgb_to_yCbCr(uint8_t *r_to_y, uint8_t *g_to_Cb, uint8_t *b_to_Cr);
void image_yCbCr_to_rgb  (uint8_t *y_to_r, uint8_t *Cb_to_g, uint8_t *Cr_to_b);

void block_downsample420(uint8_t *UINT8_BLOCK);

void block_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK);
void block_inv_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK);

void general_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_WIDTH, size_t BLOCK_HEIGHT);

void block_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK);
void block_inv_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK);

void block_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX);
void block_inv_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX);

int block_encode(OUTSTREAM* out, int16_t *INT16_SEQUENCE, int16_t *PREV_DC,
 const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS,
 const uint16_t *AC_CODES, const uint8_t *AC_VALUES, const uint8_t *AC_LENGTHS);

uint8_t min_bits(uint16_t n);
int min_bits_abs(int16_t n);
int min_bits_code(uint16_t n);

// coders
bool search_codes(INSTREAM *in, uint8_t *rrrrssss, const uint16_t *CODES, const uint8_t *VALUES, const uint8_t *LENGTHS, size_t CODES_NUMBER);
bool get_code(uint16_t *code, uint8_t rrrrssss, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER);
void encode_to_cache(uint16_t CODE, uint32_t *curr_cache, uint32_t *next_cache, int *dtr, bool min_two);

void DATA_PACKET_pack(DATA_PACKET *data, int16_t VAL, uint8_t zeros);
bool DATA_PACKET_encode(DATA_PACKET *data, const uint16_t *CODES, const uint8_t *VALUES, const uint8_t *LENGTHS, size_t N_CODES);

RGB_IMAGE *new_RGB_IMAGE(uint16_t width, uint16_t height);
void delete_RGB_IMAGE(RGB_IMAGE *img);



int OUTSTREAM_push(OUTSTREAM *out, uint16_t data, int bits);
void OUTSTREAM_reset_bytes(OUTSTREAM *out);
OUTSTREAM *new_OUTSTREAM(const char* filename, int buffer_bytes);
int delete_OUTSTREAM(OUTSTREAM *out);

int INSTREAM_pull_1bit(INSTREAM *in);
int INSTREAM_pull(INSTREAM *in, uint16_t *data, int bits);
void INSTREAM_reset_bytes(INSTREAM *in);
INSTREAM *new_INSTREAM(const char* filename, int buffer_bytes);
int delete_INSTREAM(INSTREAM *in);

int block_decode(INSTREAM* in, int16_t *INT16_SEQUENCE, int16_t* PREV_DC,
 const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS,
 const uint16_t *AC_CODES, const uint8_t *AC_VALUES, const uint8_t *AC_LENGTHS);
void decode_data(INSTREAM *in, uint8_t rrrrssss, int *ssss, int *rrrr, uint16_t *val);
bool write_data(int16_t *INT16_SEQUENCE, bool is_dc, int idx, int ssss, int rrrr, uint16_t val);

void print_32bits(uint32_t data);
void print_16bits(int16_t data);


int decode_image(char* filename, RGB_IMAGE *image);
int decode_slice(INSTREAM *in, uint8_t *slice, uint16_t WIDTH, uint16_t HEIGHT, bool is_luminance);
int encode_image(char* filename, RGB_IMAGE *image);
int encode_slice(OUTSTREAM *out, uint8_t *slice, uint16_t WIDTH, uint16_t HEIGHT, bool is_luminance);


void print_ubits(uint8_t n);
#endif