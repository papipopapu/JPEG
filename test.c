
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <math.h>
uint8_t min_bits(uint16_t n) {
    // min bits to hold the abolute value of int16 //
    int i;
    if (n < 0) n = -n;
    uint8_t count = 0;
    for (i = 15; i >= 0; i--) {
        if ((n >> i) & 1) return 16 - count;
        count++;
    }
    return 0;
}

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
typedef struct DATA_NODE {
    /* Temporary struct to pack both AC and DC components */
    uint8_t rrrrssss; // packed here
    int16_t VAL; // prob less than 16 bits, min bits will be packed and then recasted to an int16 to be interpreted
    struct DATA_NODE *next; // next pack
} DATA_NODE;
DATA_NODE *new_DATA_NODE()
{
    DATA_NODE* fetus = NULL;
    fetus = (DATA_NODE*)malloc(sizeof(DATA_NODE));
    fetus -> next = NULL;
    return fetus;
}

void free_DATA_NODE_list(DATA_NODE* head)
{
    DATA_NODE* temp;
    while (head != NULL)
    {
        temp = head;
        head = head -> next;
        free(temp);
    }
}
void connect_DATA_NODE(DATA_NODE **prev, DATA_NODE **next, DATA_NODE **head) {
    if (*prev) {(*prev) -> next = *next;} // if *prev is not null, then proceed as usual
    else {free(*head); *head = *next;} // if prev is null, then that means next is the head
    *prev = *next;
}


void pack_DATA_NODE(DATA_NODE *node, int8_t zeros, int16_t VAL)
{ 
    uint8_t isNeg = 0, minBits = min_bits_abs(VAL);
    if (VAL < 0) {VAL = -VAL;  isNeg = 1;}
    node -> rrrrssss = zeros; // number of previous zeros
    node -> rrrrssss <<= 4; // shift to the left 4 bits
    node -> rrrrssss |= minBits; // add on the other side the minimum bits to represent VAL - 1 (even though we know
    // we have eliminated the first bit, so bits = 1 -> 1, bits = 0 -> 0))
    node -> VAL = VAL; // add value on the right of VAL 
    node -> VAL <<= (15 - minBits) + 1; // shift VAL to the left as far as we can, leaving one zero for sign (the +1 is since we dont need leading 1)
    //This means that the maximum value of VAL is +-32767 !! TOTAL_BITSIZE += min_bits_abs(val);
    node -> VAL |= (isNeg << 15); // add the sign
}
bool search_codes(uint16_t compare_base, const uint16_t *CODES, int *bits_read, uint8_t *MATCHES, size_t CODES_NUMBER) {
    /* Checks if there are any matches of any amount of crecent digits of the compare base inside the CODES provided. */
    int k; *bits_read = 1; // we are goind to read at least 2 digits
    uint16_t compare_deriv;
    while ((*bits_read) <= 16) { // while we are not at the end of the group in focus and we haven't found the code
        (*bits_read)++; 
        compare_deriv = compare_base >> (16 - *bits_read); // read progressibely more digits of the code
        for (k=0; k < CODES_NUMBER; k++) { // check if code exists in our list
            if (compare_deriv == CODES[k]) {
                MATCHES[k]++;
                return true;
            }
        }  
    } return false;
}
int DECODE_BIN_16_64(FILE *file, const uint16_t *CODES, uint8_t *MATCHES, size_t CODES_NUMBER, size_t CACHES_TO_READ) {
    /* 
    Reads binary file with binary data, and compares it with a list of prefix binary codes (no code is the prefix of another),
    starting from size of 2 bits. That is, the code reads a minimum of 2 bits, but the integers correspoding to the codes 
    can have 1 bit (integers 0 and 1), as long as they are have been encoded like "00" and "01" respectively.
    
    Args:
    * file: file pointer to the binary file in "wb" mode
    * CODES: array of ints correponding decimal interpretation of the binary codes
    * MATCHES: output array of ints correponding to the number of times the code with the same index has been found in the file
    * CODES_NUMBER: number of elements in CODES and MATCHES
    * CACHES_TO_READ: number of 64 bit chunks to read from the file
    * 
    Returns:
    * 0 if no error occurred, 1 if there is not enough data in the file to read anything, -1 if some non-matching sequence was found.
    */
    
//  
//                 ______________   
//    ____________|_____________|______________   Caches are made of 64 bit integers, and a focus_block of 16 bits
//   |            |      |      |             |   over them, reading their bites, accounting for when the focus_block 
//   | curr_cache |      |      | next_cache  |   is contained fully inside curr_cache, and when it as slid between the
//   |____________|______|______|_____________|   two caches. The code is executed until a number of caches has been read,
//                | focus_block |                 or an error occured.
//                |_____________|         
//               
//   |____________|
//        dtr     
//
//  dtr is the displacement of the right of the focus_block from the left side of theL current_cache, 
//  bits is the bit size of the last code read, els is the elements read from the file (whenever we do so),
//  and caches_read is the caches read :D
    int dtr = 0, bits = 69, caches_read = 0, els;

//  create caches
    uint64_t curr_cache, next_cache; 
    uint16_t focus_block, deriv_block; 

//  is the focus_block fully contained inside the current_cache (= slipping)?
    bool SLIP = false;

//  if no block can be read, throw error 1
    els = fread(&curr_cache, 8, 1, file); 
    if (els == 0) return 1;

//  while we haven't read enough caches
    while (caches_read < CACHES_TO_READ) { 

//      slipping if it is displaced to the righ further than "sizeof(cache)=64 bits - sizeof(focus_block)=16 bits
        SLIP =  dtr > 48; 

//      if focus_block is slipping, special course of action
        if (SLIP) {

//          try to read next cache into next_cache (smart name)
            els = fread(&next_cache, 8, 1, file);

//          if no further caches can be read, then we set the next_cache to 0, and that works for us
            if (els == 0) next_cache = 0;
            while(SLIP) {             

//              get the focus_block combining both caches    2x64 -16 - dtr                  dtr + 64 - 16
                focus_block = 0; focus_block |= (next_cache >> (112 - dtr)) | (curr_cache << (dtr - 48));

//              search for matches inside the provided codes, return -1 if no match is found, add the bits of 
//              code read to the displacement
                if (!search_codes(focus_block, CODES, &bits, MATCHES, CODES_NUMBER)) return -1; dtr += bits;

//              now slip is true while we are slipping, and it stops when the whole focus block slid through to 
//              the next cache
                SLIP = dtr < 64;

//          when the block has slid through, we either return 0 if it was the last cache, or we go on and take
//          the next_cache as our new curr_cache
            } if (els == 0) return 0; curr_cache = next_cache; caches_read++; dtr = 0;

//      if focus_block is not slipping, we  proceed as usual and read the focus block entirely from curr_cache
        } else {                
            
//          get focus_block by trimming curr_cache     64 - 16 -dtr
            focus_block = 0; focus_block |= curr_cache >> (48 - dtr);

//          search for matches inside the provided codes, return -1 if no match is found, add the bits of 
//          code read to the displacement
            if (!search_codes(focus_block, CODES, &bits, MATCHES, CODES_NUMBER)) return -1; dtr += bits;   
        }
    } 

//  easy
    return 0;
}
uint16_t read_nbits(uint16_t VAL, uint8_t n) {
    /* 
    Reads n bits from a 16 bit integer, and returns the result.
    
    Args:
    * VAL: 16 bit integer to read from
    * n: number of bits to read
    
    Returns:
    * The result of reading n bits from VAL
    */
    return VAL >> (16 - n);
}
void encode_CODE(uint16_t CODE, uint64_t *curr_cache, uint64_t *next_cache, uint8_t disp) {
    *next_cache = 0; // curr_cache = [{prev_codes}0000000]
    uint64_t bigger = CODE;
    int bits = min_bits(CODE); bits += (CODE == 1 | CODE == 0) ? 1 : 0;
    int SLIP = disp + bits - 64;
    if (SLIP > 0) {
        *curr_cache |= bigger >> SLIP; // mira esto bien no esta bien .D
        *next_cache |= bigger << SLIP;
    } else {
        *curr_cache |= bigger << -SLIP;
    }
    


}
void encode_VALUE(uint16_t VALUE, uint64_t curr_cache, uint64_t next_cache) {

}
int ENCODE_DC(FILE *file, DATA_NODE *NODES, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER, size_t VALUES_NUMBER) {
    DATA_NODE *node = NODES;
    uint8_t value = node -> VAL;
    uint64_t curr_cache = 0, next_cache = 0;
    int curr_dtr = 0, prev_dtr, bits = 0;
    if (!NODES) return -1;
    while (node -> next) {
        curr_dtr += node -> rrrrssss;

        if (curr_dtr > 48) {



            fwrite(&curr_cache, 8, 1, file);
            curr_cache = next_cache;
            curr_cache -= 64; prev_dtr = curr_dtr;
        }
        
        value = node -> VAL;
        fwrite(&curr_cache, 8, 1, file);
        node = node -> next;
    } // reache last node

}
    // file in append mode

    

//            printf("DTR: %d, SLIP: %d, focus_block: ", dtr, SLIP); print_bits(focus_block); printf("\n");

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
    n = 0;
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
    int success = DECODE_BIN_16_64(file, codes, matches, 8, 2);
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