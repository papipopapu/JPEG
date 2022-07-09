
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
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
void print_bits(uint16_t n)
{
    int i;
    for (i = 15; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
}
void print_bits64(uint64_t n)
{
    int i;
    for (i = 63; i >= 0; i--) {
        printf("%ld", (n >> i) & 1);
    }
}
bool search_codes(uint16_t compare_base, const uint16_t *CODES, int *bits_read, uint8_t *MATCHES, size_t CODES_NUMBER){
    int k, bits = 1; // we are goind to read at least 2 digits
    bool code_found = false;
    uint16_t compare_deriv;
    while (bits <= 16 && !code_found) { // while we are not at the end of the group in focus and we haven't found the code
        bits++; 
        compare_deriv = compare_base >> (16 - bits); // read progressibely more digits of the code
        for (k=0; k<CODES_NUMBER; k++) { // check if code exists in our list
            if (compare_deriv == CODES[k]) {
                code_found = true;
                MATCHES[k]++;
            }
        }  
    }
    *bits_read = bits;
    return code_found;
}
int DECODE_BIN_16_64(FILE *file, const uint16_t *CODES, uint8_t *MATCHES, size_t CODES_NUMBER, size_t CACHES_TO_READ) {
    int dtr = 0, els, bits = 69, caches_read = 0;
    uint64_t curr_cache, next_cache; 
    uint16_t focus_block, deriv_block; 
    bool SLIP = false;// util, number of bits inside last code read, codes read, wether matching code has been found
    els = fread(&curr_cache, 8, 1, file); // initialize initial left and right groups
    if (CACHES_TO_READ == 0) return 1;
    if (els == 0) return 2;
    while (caches_read < CACHES_TO_READ) { 
        SLIP =  dtr > 48; // (dtr + 16) > 64 
        if (SLIP) {
            els = fread(&next_cache, 8, 1, file);
            if (els == 0) next_cache = 0;
            while(SLIP) {                                    // 2x64 -16 -dtr                // dtr + 64 - 16
                focus_block = 0; focus_block |= (next_cache >> (112 - dtr)) | (curr_cache << (dtr - 48));
                if (!search_codes(focus_block, CODES, &bits, MATCHES, CODES_NUMBER)) return -1; dtr += bits;
                SLIP = dtr < 64;
            } if (els == 0) return 0; curr_cache = next_cache; caches_read++; dtr = 0;
        } else {                                        // 64 - 16 -dtr
            focus_block = 0; focus_block |= curr_cache >> (48 - dtr);
            printf("DTR: %d, SLIP: %d, focus_block: ", dtr, SLIP); print_bits(focus_block); printf("\n");
            if (!search_codes(focus_block, CODES, &bits, MATCHES, CODES_NUMBER)) return -1; dtr += bits;   
        }
    } 
    return 0;

}
int decode_uniquecodes_bin(FILE *file, const uint16_t *CODES, uint8_t *MATCHES, size_t CODE_NUMBER, size_t CODES_TO_READ, size_t CACHE_BLOCKS_BYTES) {
 if (CODES_TO_READ == 0) return -1;
    // group = 16 bit number fromed by taking a 16 bits of the binary file 
    uint64_t left, right;
    uint16_t compare_base, compare_deriv;// left group, right group, group in focus, util
    int dtr = 0, els, CACHE_BLOCKS_BITS = CACHE_BLOCKS_BYTES*8;//  displacement of compare_base to the right with respect to the left most side of left group
    int k, bits = 0, codes_read = 0;
    int bits_left = 1; 
    bool code_found = false, end_reached = false;// util, number of bits inside last code read, codes read, wether matching code has been found
    els = fread(&left, CACHE_BLOCKS_BYTES, 1, file); // initialize initial left and right groups
    if (els == 0) return -1;
    els = fread(&right, CACHE_BLOCKS_BYTES, 1, file);
    if (els == 0) right = 0;
    do { // only read CODES_TO_READ codes
        if (dtr  >= (CACHE_BLOCKS_BITS-16)) { // if the compare_base is completely out of the left group
            left = right; // move the left group to the right group
            els = fread(&right, CACHE_BLOCKS_BYTES, 1, file);
            if (els == 0) {
                end_reached = true;
                right = 0;
                bits_left = CACHE_BLOCKS_BITS - dtr;
            }
            dtr -= CACHE_BLOCKS_BITS; // reset the displacement taking into account the slip into the new right group
            bits = dtr; // i0 = displacement to the right after the reset
        }
        compare_base = 0;// 
        compare_base |= (dtr > (CACHE_BLOCKS_BITS - 16)) ? (right >> (2*CACHE_BLOCKS_BITS - 16 - dtr)) | (left << (16 - CACHE_BLOCKS_BITS + dtr)) : (left >> (CACHE_BLOCKS_BITS - 16 - dtr)); // according the the amount of bits inside the las code read, get next group in focus
        bits = 1; // we are goind to read at least 2 digits
        while (bits <= 16 && !code_found) { // while we are not at the end of the group in focus and we haven't found the code
            bits++; // see pre-previous comment for explanation on order of execution
            compare_deriv = compare_base >> (16 - bits); // read progressibely more digits of the code
            for (k=0; k<CODE_NUMBER; k++) { // check if code exists in our list
                if (compare_deriv == CODES[k]) {
                    code_found = true;
                    MATCHES[k]++;
                }
            }  
        }
        if (!code_found) {return 1;} // no corresponding code found, something went wrong...
        codes_read++; code_found = false; dtr += bits; if (end_reached) bits_left -= bits; // updates
    } while (codes_read < CODES_TO_READ &&  bits_left > 0);
    return 0;
}
const uint16_t codes[] = {
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
    uint8_t matches[] = {0,0,0,0,0,0,0,0};
    size_t len = 0;
    FILE *file = fopen("data.bin", "wb");
    uint64_t n = 15271705341792288768, n2 = 0; // 26 and 1, 11010 01 0 0000 0000
    len = fwrite(&n, sizeof(uint64_t), 1, file);
    if (len != 1) {
        printf("error writing to file\n");
        return 1;
    }
    fclose(file);

    file = fopen("data.bin", "rb");
    fread(&n2, sizeof(uint64_t), 1, file);
    printf("Truth: %lu\n", n2);
    printf("true bits: "); print_bits64(n2); printf("\n");
    fclose(file);

    uint32_t n3 = 0;
    file = fopen("data.bin", "rb");
    fread(&n3, sizeof(uint32_t), 1, file);
    printf("TEST: %u\n", n3);
    printf("TESTA bits: "); print_bits64(n3); printf("\n");
    fclose(file);









    file = fopen("data.bin", "rb");
    int success = DECODE_BIN_16_64(file, codes, matches, 8, 1);
    fclose(file);

    printf("Success?: %d\n", success);

    for (int i = 0; i < 8; i++) {
        printf("%d: %d\n", codes[i], matches[i]);
    }



    

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

    uint64_t n = 257;


    
    if (!file) return 1;
    els = fwrite(&n, sizeof(uint64_t), 1, file);
    fclose(file);
    if (els == 0) return 2;

    // read
    file = fopen("data.bin", "rb");
    }
    uint8_t n2;
    
    if (!file) return 1;
    els = fread(&n2, sizeof(uint8_t), 1, file);
    fclose(file);
    if (els == 0) return 2;

    printf("%d\n", n2);
    printf("%d\n", 1 << 3);
    
    */