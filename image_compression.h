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
extern const u_int8_t ZIGZAG_IDX_8_8[64];

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
void connect_DATA_NODE(DATA_NODE **prev, DATA_NODE **next, DATA_NODE **head);

// block_process
void get_block(u_int8_t *IMAGE, u_int8_t *UINT8_BLOCK, size_t BLOCK_SIZE, size_t IMG_WIDTH, size_t IMG_HEIGHT, size_t I0, size_t J0);

void block_rgb_to_yCbCr(u_int8_t *r_to_y, u_int8_t *g_to_Cb, u_int8_t *b_to_Cr, size_t BLOCK_DIM);
void block_yCbCr_to_rgb  (u_int8_t *y_to_r, u_int8_t *Cb_to_g, u_int8_t *Cr_to_b, size_t BLOCK_DIM);

void block_downsample420(u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM);

void block_dct(u_int8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_DIM);
void block_inv_dct(u_int8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_DIM);

void general_dct(u_int8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_WIDTH, size_t BLOCK_HEIGHT);

void block_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_DIM);
void block_inv_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK,  size_t BLOCK_DIM);

void block_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, size_t BLOCK_DIM, const u_int8_t *SERIAL_IDX);
void block_inv_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, size_t BLOCK_DIM, const u_int8_t *SERIAL_IDX);

u_int8_t min_bits_abs(int16_t n);
void blocks_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, size_t BLOCK_DIM, size_t BLOCK_NUMBER);
void block_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, size_t BLOCK_DIM, bool IS_FIRST);
void block_process_one(bool isY, u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM, DATA_NODE **AC_HEAD, DATA_NODE **DC_HEAD);

void blocks_encode(FILE* file, DATA_NODE *AC_DATA_NODES, DATA_NODE *DC_DATA_NODES, int16_t IMG_WIDTH, int16_t IMG_HEIGHT);

// {uint8_block} -> dct -> {float_block} -> quantize -> {int8_block} -> serialize_reorder -> {int8_sequence} -> huffman -> {huffman}
// 16 one bits + img width(16 bits) + img height(16 bits) + 16 one bits + dc codes
#endif