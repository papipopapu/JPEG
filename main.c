#include <stdio.h>
#include <stdlib.h>
#include "image_compression.h"
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

    u_int8_t max_block[64];
    max_block[0] = 0;
    for (int k=1; k<64; k++) {
        max_block[k] = max_block[k-1] + 1;
    }

    int16_t h[64];
    u_int8_t aprox_block[n*m];
    get_block(max_block, aprox_block, 8, 8, 8, 0, 0);


    DATA_NODE* AC = new_DATA_NODE();  DATA_NODE* DC = new_DATA_NODE();
    block_process_one(true, aprox_block, 8, &AC, &DC);
    free_DATA_NODE_list(AC); free_DATA_NODE_list(DC);

    /*
    this seg-faults:
        DATA_NODE** AC;  DATA_NODE** DC;
        block_process_one(true, aprox_block, 8, AC, DC);
        free_DATA_NODE_list(*AC); free_DATA_NODE_list(*DC);
    */



    return 0;
}