#include "image_compression.h"


const uint16_t codes[] = {
0, // 00 special cases
1, // 01%d\n", 0b1101010111101101, 
1015,
4086,
32706,
65420,
};




int main () {
    
    //printf("Test basic use, writing [1101010], [1111], [0], [110101001] into a cache of 3 bytes.\n Then reading 16 bytes as [1101010111101101].\n");

    OUTSTREAM *out = new_OUTSTREAM("test.bin", 1);  
    OUTSTREAM_push(out, 0b110000011, 9);             
    delete_OUTSTREAM(out);

    FILE *file = fopen("test.bin", "rb");
    uint8_t *truth = malloc(sizeof(uint8_t)*2);
    fread(truth, 1, 2, file);
    for (int i = 0; i < 2; i++) {
        print_ubits(truth[i]);
    }
    printf("\n");
    fclose(file);
    free(truth);

    
    uint16_t num = 0;
    INSTREAM *in = new_INSTREAM("test.bin", 1);
    if (!in) {
        printf("Error: could not open file.\n");
        return 1;
    }
    INSTREAM_pull(in, &num, 9);
    delete_INSTREAM(in);


    printf("Read: %d\n", num);



 




    


    return 0;
}





