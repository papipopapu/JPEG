
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

/* for AC codes bits = 8 bits for DC codes bits = 4 bits 
1- assign weights to DATA_NODEs and store them in HUFFMN_NODEs
2- build the HUFFMAN tree
3- traverse the HUFFMAN and create table of codes
4- encode the DATA_NODEs with the table of codes
*/
int from_binary_to_decimal(char* binary_string) {
    return (int) strtol(binary_string, NULL, 2);
}
void print_bits(int16_t n)
{
    int i;
    for (i = 15; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
}

int decode_uniquecodes_bin(FILE *file, const u_int16_t *CODES, u_int8_t *MATCHES, size_t CODE_NUMBER, size_t CODES_TO_READ) {
    u_int16_t left, right, compare_base, compare_deriv;
    size_t k, i = 0, codes_read = 0; bool code_found = false;
    size_t displace_to_right = 0;
    fread(&left, sizeof(u_int16_t), 1, file);
    fread(&right, sizeof(u_int16_t), 1, file);
    while (codes_read < CODES_TO_READ) {
        if (displace_to_right  % 16 == 0 && displace_to_right != 0) {
            left = right;
            fread(&right, sizeof(u_int16_t), 1, file);
        }
        compare_base = (right >> 16 - i) | (left << i);
        i = 1; // we are goind to read at least 2 digits
        while (i <= 16 && !code_found) {
            i++;
            compare_deriv = compare_base >> 16 - i;
            for (k=0; k<CODE_NUMBER; k++) {
                if (right == CODES[k]) {
                    code_found = true;
                    MATCHES[k]++;
                }
            }  
        }
        if (!code_found) {return 1;} // no corresponding code found, something went wrong...
        codes_read++; code_found = false; displace_to_right += i;
    }
    return 0;

}

const u_int16_t codes[] = {
0, // 00 special cases
1, // 01
26,
247,
1015,
4086,
32706,
65420,
};

int main () {

    
    FILE *file = fopen("data.bin", "wb");
    fwrite(&n, sizeof(u_int16_t), 1, file);

    u_int16_t matches[] = {0, 0, 0, 0, 0, 0, 0};
    
    decode_uniquecodes_bin("data.bin", codes, matches, sizeof(codes)/sizeof(u_int16_t), 2);

    //print_bits(n);

    return 0;
}







/*
    FILE *file;
    FILE *file_out;

    
    
    size_t els;

    // write
    if (true){
    file = fopen("data.bin", "wb");

    u_int64_t n = 257;


    
    if (!file) return 1;
    els = fwrite(&n, sizeof(u_int64_t), 1, file);
    fclose(file);
    if (els == 0) return 2;

    // read
    file = fopen("data.bin", "rb");
    }
    u_int8_t n2;
    
    if (!file) return 1;
    els = fread(&n2, sizeof(u_int8_t), 1, file);
    fclose(file);
    if (els == 0) return 2;

    printf("%d\n", n2);
    printf("%d\n", 1 << 3);
    
    */