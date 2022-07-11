#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
bool is_little_endian() {
    uint16_t a = 1;
    return *(char *)&a == 1;
}
typedef struct OUTSTREAM {
    const char* filename;
    FILE *file;
    int written_bits, written_bytes, buffer_bytes;
    char *buffer;
    bool eof;
} OUTSTREAM;

 pushto_OUTSTREAM(OUTSTREAM *out, uint data, int bits) {
    // bits are the number of bits to encode from least to most significant
    if (bits < 1 || bits > 32) return 1;
    if (out->eof) return -1;
    data &= (bits == 32) ? ~0 : (1 << bits) - 1; // remove trash data, but if whole int is needed, we cant bitshift by 0
    int slip = bits + out->written_bits - 8;
    
    if (slip > 0) {
        out->buffer[out->written_bytes] |= data >> +slip;
        if (out->written_bytes == (out->buffer_bytes-1)) { // there has been slip, and at we had alredy written up to the second to last byte
            fwrite(out->buffer, 1, out->buffer_bytes, out->file);
            out->written_bytes = -1;    
        }   out->written_bits = 0;
            out->written_bytes++;

        data &= (1 << slip) - 1; // remove alredy coded part
        pushto_OUTSTREAM(out, data, slip);
    } else if (slip < 0) {
        out->buffer[out->written_bytes] |= data << -slip;
        out->written_bits += bits; 
    } else {
        out->buffer[out->written_bytes] |= data;
        out->written_bits = 0; out->written_bytes++;
    } 
    printf("Cache: \n");
    for (int i = 0; i < out->buffer_bytes; i++) {
        print_ubits(out->buffer[i]);
    }
    printf("\n");
    return 0;
}
OUTSTREAM *new_OUTSTREAM(const char* filename, int buffer_bytes) {
    if (buffer_bytes < 1) return NULL;
    OUTSTREAM *out = (OUTSTREAM *)malloc(sizeof(OUTSTREAM));
    out->filename = filename;
    out->file = fopen(filename, "wb");
    out->written_bits = 0;
    out->written_bytes = 0;
    out->eof = false;
    out->buffer_bytes = buffer_bytes;
    out->buffer = (char *)calloc(buffer_bytes, 1); 
    //if (fread(out->buffer, 1, buffer_bytes, out->file)==0) {fclose(out->file); free(out); return NULL;}
    return out;
}

int delete_OUTSTREAM(OUTSTREAM *out) {
    if (out == NULL) return -1;
    if (out->written_bits != 0 | out->written_bytes!=0) fwrite(out->buffer, 1, out->buffer_bytes, out->file); // flush buffer if something has been written
    fclose(out->file);
    free(out->buffer);
    free(out);
    return 0;
}

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





