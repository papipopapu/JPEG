#include <stdio.h>
#include <stdlib.h>
#include "image_compression.h"

void print_RGB_IMAGE(RGB_IMAGE *img) {
    printf("R:\n");
    for (int i = 0; i < img->HEIGHT; i++) {
        for (int j = 0; j < img->WIDTH; j++) {
            printf("%d, ", img->r[i*img->WIDTH+j]);
        }
        printf("\n");
    }
    printf("G:\n");
    for (int i = 0; i < img->HEIGHT; i++) {
        for (int j = 0; j < img->WIDTH; j++) {
            printf("%d, ", img->g[i*img->WIDTH+j]);
        }
        printf("\n");
    }
    printf("B:\n");
    for (int i = 0; i < img->HEIGHT; i++) {
        for (int j = 0; j < img->WIDTH; j++) {
            printf("%d, ", img->b[i*img->WIDTH+j]);
        }
        printf("\n");
    }
}
int main() {
    
    size_t n = 8, m = 8;
    int16_t prev_dc = 0;
    uint16_t width, height;
    int16_t import_block[64];
    uint8_t uint8_block[64] = {
        52, 55, 61, 66, 70, 61, 64, 73,
        63, 59, 55, 90, 109, 85, 69, 72, 
        62, 59, 68, 113, 144, 104, 66, 73,
        63, 58, 71, 122, 154, 106, 70, 69,
        67, 61, 68, 104, 126, 88, 68, 70,
        79, 65, 60, 70, 77, 68, 58, 75, 
        85, 71, 64, 59, 55, 61, 65, 83,
        87, 79, 69, 68, 65, 76, 78, 94
    };
    RGB_IMAGE *img = new_RGB_IMAGE(8, 8); 
    RGB_IMAGE *img_out = new_RGB_IMAGE(8, 8);
    memcpy(img->r, uint8_block, sizeof(uint8_t) * 64);
    memcpy(img->g, uint8_block, sizeof(uint8_t) * 64);
    memcpy(img->b, uint8_block, sizeof(uint8_t) * 64);

    encode_image("test.bin", img);
    decode_image("test.bin", img_out);

    print_RGB_IMAGE(img);
    print_RGB_IMAGE(img_out);

    // INSTREAM* in = new_INSTREAM("test.bin", 8);
    // INSTREAM_pull(in, &width, 16); INSTREAM_pull(in, &height, 16);
    // printf("Width: %d, Height: %d\n", width, height);


    // printf("Result: %d\n", block_decode(in, import_block, &prev_dc, DC_LUMINANCE_CODES, DC_VALUES, DC_LUMINANCE_LENGTHS, AC_LUMINANCE_CODES, AC_VALUES, AC_LUMINANCE_LENGTHS));
    // print_matrix(import_block);
    // prev_dc = 0;
    // printf("Result: %d\n", block_decode(in, import_block, &prev_dc, DC_LUMINANCE_CODES, DC_VALUES, DC_LUMINANCE_LENGTHS, AC_LUMINANCE_CODES, AC_VALUES, AC_LUMINANCE_LENGTHS));
    // print_matrix(import_block);
    // prev_dc = 0;
    // printf("Result: %d\n", block_decode(in, import_block, &prev_dc, DC_LUMINANCE_CODES, DC_VALUES, DC_LUMINANCE_LENGTHS, AC_LUMINANCE_CODES, AC_VALUES, AC_LUMINANCE_LENGTHS));
    // print_matrix(import_block);

    // delete_INSTREAM(in);
    delete_RGB_IMAGE(img); 
    delete_RGB_IMAGE(img_out);




    /*
    int16_t sequence[64];
    float float_block[64];
    int16_t int16_block[64];
    int16_t import_block[64];

    block_dct(uint8_block, float_block);
    block_quantize(LUMINANCE_QUANT , int16_block, float_block );
    block_serialize(int16_block, sequence,  ZIGZAG_IDX);
    

    OUTSTREAM *out = new_OUTSTREAM("out.txt", 8);
    printf("Result: %d\n", block_encode(out, sequence, 0, DC_LUMINANCE_CODES, DC_VALUES, DC_LUMINANCE_LENGTHS, AC_LUMINANCE_CODES, AC_VALUES, AC_LUMINANCE_LENGTHS));
    delete_OUTSTREAM(out);

    INSTREAM* in = new_INSTREAM("out.txt", 8);
    printf("Result: %d\n", block_decode(in, import_block, &prev_dc, DC_LUMINANCE_CODES, DC_VALUES, DC_LUMINANCE_LENGTHS, AC_LUMINANCE_CODES, AC_VALUES, AC_LUMINANCE_LENGTHS));
    delete_INSTREAM(in);

    print_matrix(import_block);
    printf("\n");
    print_matrix(sequence);
    */

    
    //printf("VAL_bits: "); print_ubits(test.VAL_bits); printf("\n");




    /*
    this seg-faults:
        DATA_NODE** AC;  DATA_NODE** DC;
        block_process_one(true, aprox_block, 8, AC, DC);
        free_DATA_NODE_list(*AC); free_DATA_NODE_list(*DC);
    */


    
    return 0;
}