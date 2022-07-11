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
// Constants
// BLOCK DIMENSIONS = 8X8
extern const float LUMINANCE_QUANT[64];
extern const float CHROMINANCE_QUANT[64];
extern const uint8_t ZIGZAG_IDX[64];
extern const uint16_t AC_LUMINANCE_HCODES[162];
extern const uint16_t AC_CHROMINANCE_CODES[162];
extern const uint8_t AC_VALUES[162];
extern const uint16_t DC_LUMINANCE_CODES[12];
extern const uint16_t DC_CHROMINANCE_CODES[12];
extern const uint8_t DC_VALUES[12];

// data_node
typedef struct DATA_NODE {
    /* Temporary struct to pack both AC and DC components */
    uint8_t rrrrssss; // packed here
    uint16_t VAL; // prob less than 16 bits, min bits will be packed and then recasted to an int16 to be interpreted
    struct DATA_NODE *next; // next pack
} DATA_NODE;
typedef struct OUTSTREAM {
    const char* filename;
    FILE *file;
    int d_bits, d_bytes, buffer_bytes;
    char *buffer;
    bool eof;
} OUTSTREAM;

typedef struct DECODER {
    const char* filename;
    FILE *file;
    int dtr;
    uint32_t curr_cache, next_cache;
    bool end_of_file;
} DECODER;

DATA_NODE *new_DATA_NODE();
void free_DATA_NODE_list(DATA_NODE* head);
void pack_DATA_NODE(DATA_NODE *node, int8_t zeros, int16_t VAL);
void connect_DATA_NODE(DATA_NODE **prev, DATA_NODE **next, DATA_NODE **head);

// block_process
void get_block(uint8_t *IMAGE, uint8_t *UINT8_BLOCK, size_t IMG_WIDTH, size_t IMG_HEIGHT, size_t I0, size_t J0);

void block_rgb_to_yCbCr(uint8_t *r_to_y, uint8_t *g_to_Cb, uint8_t *b_to_Cr);
void block_yCbCr_to_rgb  (uint8_t *y_to_r, uint8_t *Cb_to_g, uint8_t *Cr_to_b);

void block_downsample420(uint8_t *UINT8_BLOCK);

void block_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK);
void block_inv_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK);

void general_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_WIDTH, size_t BLOCK_HEIGHT);

void block_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK);
void block_inv_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK);

void block_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX);
void block_inv_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX);

uint8_t min_bits(uint16_t n);
void blocks_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, size_t BLOCK_NUMBER);
void block_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, bool IS_FIRST);
void block_process_one(bool isY, uint8_t *UINT8_BLOCK, DATA_NODE **AC_HEAD, DATA_NODE **DC_HEAD);

void blocks_encode(FILE* file, DATA_NODE *AC_DATA_NODES, DATA_NODE *DC_DATA_NODES,
    const uint16_t *DC_NEWCODES, const uint16_t *DC_OLDCODES, const uint16_t *AC_NEWCODES, const uint16_t *AC_OLDCODES,
    size_t BLOCK_NUMBER, size_t CODES_NUMBER);

// coders
bool search_codes(uint16_t compare_base, const uint16_t *CODES, int *bits_read, uint8_t *MATCHES, size_t CODES_NUMBER);
bool get_code(uint16_t *code, uint8_t rrrrssss, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER);
void encode_to_cache(uint16_t CODE, uint32_t *curr_cache, uint32_t *next_cache, int *dtr, bool min_two);
uint16_t get_bits_at(uint32_t cache, int dtr, int bits);
int ENCODE_DATA(const char* filename, DATA_NODE *NODES, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER);

uint16_t pullfrom_DECODER(DECODER *decoder, int bits);
DECODER *new_DECODER(const char* filename);
void free_DECODER(DECODER *decoder);

void pushto_ENCODER(ENCODER *encoder, uint16_t CODE, bool min_two);
ENCODER *new_ENCODER(const char* filename);
void free_ENCODER(ENCODER *encoder);





#endif