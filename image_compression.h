#ifndef IMAGE_COMPRESSION_H
#define IMAGE_COMPRESSION_H

// Includes
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>


// Constants
extern const float LUMINANCE_QUANT_MATRIX_8_8[64];
extern const float CHROMINANCE_QUANT_MATRIX_8_8[64];

// data_node
typedef struct DATA_NODE {
    /* Temporary struct to pack both AC and DC components */
    u_int8_t zeros_bitsVAL; // packed here
    int16_t VAL; // prob less than 16 bits, min bits will be packed and then recasted to an int16 to be interpreted
    struct DATA_NODE *next; // next pack
} DATA_NODE;

DATA_NODE *new_DATA_NODE();
void free_DATA_NODE_list(DATA_NODE* head);
void pack_DATA_NODE(DATA_NODE *node, int8_t zeros, int16_t VAL);

// block_process
void get_block(u_int8_t *IMAGE, u_int8_t *UINT8_BLOCK, size_t BLOCK_SIZE, size_t IMG_WIDTH, size_t IMG_HEIGHT, size_t I0, size_t J0);
void bloc_rgb_to_yCbCr(u_int8_t *r_to_y, u_int8_t *g_to_Cb, u_int8_t *b_to_Cr, size_t BLOCK_DIM);
void block_downsample420(u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM);
void block_dct(u_int8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_DIM);
void block_idct(float *FLOAT_BLOCK, u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM);
void general_dct(u_int8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_WIDTH, size_t BLOCK_HEIGHT);
void block_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_DIM);
void block_zigzag(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, size_t BLOCK_DIM);
u_int8_t min_bits_abs(int16_t n);
void block_pack(int16_t *INT16_SEQUENCE, DATA_NODE *PACKED_BLOCK_HEAD, size_t BLOCK_DIM); 
void block_process_one(bool isY, u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM, DATA_NODE *PACKED_BLOCK_HEAD);

// {uint8_block} -> dct -> {float_block} -> quantize -> {int8_block} -> zigzag_reorder -> {int8_sequence} -> huffman -> {huffman}

#endif