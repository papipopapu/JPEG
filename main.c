#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "core/image_compression.h"

void print_block(uint8_t *block) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%u ", block[i*8+j]);
        }
        printf("\n");
    }
}

int main(void) {
    uint16_t w, h;
    uint8_t *rin, *gin, *bin;
    uint8_t *rout = malloc(sizeof(uint8_t) * 64);
    uint8_t *gout = malloc(sizeof(uint8_t) * 64);
    uint8_t *bout = malloc(sizeof(uint8_t) * 64);



    for (int i = 0; i < 64; i++) {
        rout[i] = rand();
    }

    for (int i = 0; i < 64; i++) {
        gout[i] = rand();
    }
    // fill bin array with random values
    for (int i = 0; i < 64; i++) {
        bout[i] = rand();
    }

    printf("R:\n");
    print_block(rout); 
    printf("G:\n");
    print_block(gout);
    printf("B:\n");
    print_block(bout);


    encode_image("img.bin", rout, gout, bout, 8, 8);
    decode_image("img.bin", &rin, &gin, &bin, &w, &h);

    printf("R:\n");
    print_block(rin); 
    printf("G:\n");
    print_block(gin);
    printf("B:\n");
    print_block(bin);

    free(rin); free(gin); free(bin); free(rout); free(gout); free(bout);

    return 0;
}