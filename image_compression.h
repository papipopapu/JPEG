#ifndef IMAGE_COMPRESSION_H
#define IMAGE_COMPRESSION_H

// Includes
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>


// Constants
// BLOCK DIMENSIONS = 8X8
extern const float LUMINANCE_QUANT[64];
extern const float CHROMINANCE_QUANT[64];
extern const uint8_t ZIGZAG_IDX[64];
extern const uint16_t AC_LUMINANCE_HUFF[162];
extern const uint16_t AC_CHROMINANCE_HUFF[162];
extern const uint8_t AC_HUFF_IDX[162];
extern const uint16_t DC_LUMINANCE_HUFF[12];
extern const uint16_t DC_CHROMINANCE_HUFF[12];
extern const uint8_t DC_HUFF_IDX[12];

// data_node
typedef struct DATA_NODE {
    /* Temporary struct to pack both AC and DC components */
    uint8_t zeros_bitsVAL; // packed here
    int16_t VAL; // prob less than 16 bits, min bits will be packed and then recasted to an int16 to be interpreted
    struct DATA_NODE *next; // next pack
} DATA_NODE;


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

uint8_t min_bits_abs(int16_t n);
void blocks_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, size_t BLOCK_NUMBER);
void block_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, bool IS_FIRST);
void block_process_one(bool isY, uint8_t *UINT8_BLOCK, DATA_NODE **AC_HEAD, DATA_NODE **DC_HEAD);

void blocks_encode(FILE* file, DATA_NODE *AC_DATA_NODES, DATA_NODE *DC_DATA_NODES,
    const uint16_t *DC_NEWCODES, const uint16_t *DC_OLDCODES, const uint16_t *AC_NEWCODES, const uint16_t *AC_OLDCODES,
    size_t BLOCK_NUMBER, size_t CODES_NUMBER);
// {uint8_block} -> dct -> {float_block} -> quantize -> {int8_block} -> serialize_reorder -> {int8_sequence} -> huffman -> {huffman}

// 16 one bits + img width(16 bits) + img height(16 bits) + 16 one bits 
// + dc codes/values + 16ob + ac codes/values  + 16ob (for yCbCr)
// fill rest




#endif