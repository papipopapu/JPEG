#include "image_compression.h"


const uint16_t codes[] = {
0, // 00 special cases
1, // 01%d\n", 0b1101010111101101, 
1015,
4086,
32706,
65420,
};

void print_bits_uint16(uint16_t n)
{
    int i;
    for (i = 15; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
}


int main () {
    
    //printf("Test basic use, writing [1101010], [1111], [0], [110101001] into a cache of 3 bytes.\n Then reading 16 bytes as [1101010111101101].\n");

    OUTSTREAM *out = new_OUTSTREAM("test.bin", 2);  
    pushto_OUTSTREAM(out, 0b000000001111000011111111, 24);
                       
    delete_OUTSTREAM(out);

    FILE * file = fopen("test.bin", "rb");
    char *num = (char *)malloc(3);
    printf("Els read: %lu\n",fread(num, 1, 3, file));

    printf("Extract all: \n");
    for (int i = 0; i < 3; i++) {
        print_ubits((uint8_t)num[i]);
    }
    printf("\n");

    free(num);
    fclose(file);





    


    return 0;
}





