#include <stdio.h>
#include <stdlib.h>
#include "image_compression.h"
void print_matrix(int16_t *seq) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%d ", seq[i*8+j]);
        }
        printf("\n");
    }
}
int main() {
    
    size_t n = 8, m = 8;
    u_int8_t uint8_block[64] = {
        52, 55, 61, 66, 70, 61, 64, 73,
        63, 59, 55, 90, 109, 85, 69, 72, 
        62, 59, 68, 113, 144, 104, 66, 73,
        63, 58, 71, 122, 154, 106, 70, 69,
        67, 61, 68, 104, 126, 88, 68, 70,
        79, 65, 60, 70, 77, 68, 58, 75, 
        85, 71, 64, 59, 55, 61, 65, 83,
        87, 79, 69, 68, 65, 76, 78, 94
    };

    int16_t sequence[64];
    float float_block[64];
    int16_t int16_block[64];

    block_dct(uint8_block, float_block);
    block_quantize(LUMINANCE_QUANT , int16_block, float_block );
    block_serialize(int16_block, sequence,  ZIGZAG_IDX);
    print_matrix(sequence);

    OUTSTREAM *out = new_OUTSTREAM("out.txt", 8);
    block_encode(out, sequence, 0, DC_LUMINANCE_CODES, DC_VALUES, AC_LUMINANCE_CODES, AC_VALUES);
    delete_OUTSTREAM(out);

    
    //printf("VAL_bits: "); print_ubits(test.VAL_bits); printf("\n");




    /*
    this seg-faults:
        DATA_NODE** AC;  DATA_NODE** DC;
        block_process_one(true, aprox_block, 8, AC, DC);
        free_DATA_NODE_list(*AC); free_DATA_NODE_list(*DC);
    */



    return 0;
}