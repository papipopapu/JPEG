
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <inttypes.h>
#include <math.h>
#define min(a,b) \
({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a < _b ? _a : _b; })


uint8_t min_bits(uint16_t n) {  
    if (n == 0) return 1;
    int i;
    uint8_t count = 0;
    for (i = 15; i >= 0; i--) {
        if ((n >> i) & 1) return 16 - count;
        count++;
    }
    return 0;
}

/* for AC codes bits = 8 bits for DC codes bits = 4 bits print_bits32
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
void print_bits32(uint32_t n)
{
    int i;
    for (i = 31; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
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
    uint8_t isNeg = 0, minBits; 
    if (VAL < 0) {VAL = -VAL;  isNeg = 1;}
    minBits = (VAL == 0) ? 0 : min_bits(VAL);
    node -> rrrrssss = zeros; // number of previous zeros
    node -> rrrrssss <<= 4; // shift to the left 4 bits
    node -> rrrrssss |= minBits; // add on the other side the minimum bits to represent VAL - 1 (even though we know
    // we have eliminated the first bit, so bits = 1 -> 1, bits = 0 -> 0))
    node -> VAL = VAL; // add value on the right of VAL 
    node -> VAL &= ~(1 << (minBits-1)); // mask off the bits that are not needed
    node -> VAL |= isNeg << (minBits-1); // shift VAL to the left as far as we can, leaving one zero for sign (the +1 is since we dont need leading 1)
    //This means that the maximum value of VAL is +-32767 !! TOTAL_BITSIZE += min_bits_abs(val);
     // add the sign
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
int DECODE_BIN_16_32(FILE *file, const uint16_t *CODES, uint8_t *MATCHES, size_t CODES_NUMBER, size_t CACHES_TO_READ) {
    /* 
    Reads binary file with binary data, and compares it with a list of prefix binary codes (no code is the prefix of another),
    starting from size of 2 bits. That is, the code reads a minimum of 2 bits, but the integers correspoding to the codes 
    can have 1 bit (integers 0 and 1), as long as they are have been encoded like "00" and "01" respectively.
    
    Args:
    * file: file pointer to the binary file in "wb" mode
    * CODES: array of ints correponding decimal interpretation of the binary codes
    * MATCHES: output array of ints correponding to the number of times the code with the same index has been found in the file
    * CODES_NUMBER: number of elements in CODES and MATCHES
    * CACHES_TO_READ: number of 32 bit chunks to read from the file
    * 
    Returns:
    * 0 if no error occurred, 1 if there is not enough data in the file to read anything, -1 if some non-matching sequence was found.
    */
    
//  
//                 ______________   
//    ____________|_____________|______________   Caches are made of 32 bit integers, and a focus_block of 16 bits
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
    uint32_t curr_cache, next_cache; 
    uint16_t focus_block, deriv_block; 

//  is the focus_block fully contained inside the current_cache (= slipping)?
    bool SLIP = false;

//  if no block can be read, throw error 1
    els = fread(&curr_cache, 8, 1, file); 
    if (els == 0) return 1;

//  while we haven't read enough caches
    while (caches_read < CACHES_TO_READ) { 

//      slipping if it is displaced to the righ further than "sizeof(cache)=32 bits - sizeof(focus_block)=16 bits
        SLIP =  dtr > 16; 

//      if focus_block is slipping, special course of action
        if (SLIP) {

//          try to read next cache into next_cache (smart name)
            els = fread(&next_cache, 8, 1, file);

//          if no further caches can be read, then we set the next_cache to 0, and that works for us
            if (els == 0) next_cache = 0;
            while(SLIP) {             

//              get the focus_block combining both caches    2x32 -16 - dtr                  dtr - 32 + 16
                focus_block = 0; focus_block |= (next_cache >> (48 - dtr)) | (curr_cache << (dtr - 16));

//              search for matches inside the provided codes, return -1 if no match is found, add the bits of 
//              code read to the displacement
                if (!search_codes(focus_block, CODES, &bits, MATCHES, CODES_NUMBER)) return -1; dtr += bits;

//              now slip is true while we are slipping, and it stops when the whole focus block slid through to 
//              the next cache
                SLIP = dtr < 32;

//          when the block has slid through, we either return 0 if it was the last cache, or we go on and take
//          the next_cache as our new curr_cache
            } if (els == 0) return 0; curr_cache = next_cache; caches_read++; dtr = 0;

//      if focus_block is not slipping, we  proceed as usual and read the focus block entirely from curr_cache
        } else {                
            
//          get focus_block by trimming curr_cache     32 - 16 -dtr
            focus_block = 0; focus_block |= curr_cache >> (16 - dtr);

//          search for matches inside the provided codes, return -1 if no match is found, add the bits of 
//          code read to the displacement
            if (!search_codes(focus_block, CODES, &bits, MATCHES, CODES_NUMBER)) return -1; dtr += bits;   
        }
    } 

//  easy
    return 0;
}

bool get_code(uint16_t *code, uint8_t rrrrssss, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER) {
    int i;
    for (i=0; i < CODES_NUMBER; i++) {
        if (rrrrssss == VALUES[i]) {*code = CODES[i]; return true;}
    } return false;
}
void encode_16_32(uint16_t CODE, uint32_t *curr_cache, uint32_t *next_cache, int *dtr, bool HUFFMAN_MODE) {
    //*next_cache = [00000000...];  curr_cache = [{prev_codes}{0000000}...]
    uint32_t bigger = CODE; int bits, slip;
    bits = min_bits(CODE); if (HUFFMAN_MODE) {bits += ((CODE == 1 || CODE == 0) ? 1 : 0);}
    slip = *dtr + bits - 32;
    if (slip > 0) {
        *curr_cache |= bigger >> slip; 
        *next_cache |= bigger << (32 - slip);
    } else if (slip == 0) {
        *curr_cache |= bigger;
    } else {
        *curr_cache |= bigger << -slip;
    } *dtr += bits;
}

uint16_t get_bits_at(uint32_t cache, int dtr, int bits) {
    /*  Get a number of bits at a position dtr (displacement to the right) from a uint32_t cache.
        If bits produce overflow, they are truncated to the maximum amount of bits that fit.
        If dtr is bigger or equal to 32, dtr < 0 or the number of bits asked is 0, the function returns 0.

        Args:
        * cache: the cache to read from
        * dtr: the displacement to the right
        * bits: the number of bits to read
        Returns:
        * The bits read from the cache
    
    [0,0,0,{1,0,1,1,1},...] size = 32 bits
       |____| |________|
        dtr      bits
    */  
    if (dtr >= 32 || dtr < 0 || bits == 0) return 0;
    bits = min(bits, 32 - dtr);
    int desp; uint16_t ret;
    desp = 32 - dtr - bits;
    ret = (desp == 0) ? cache : cache >> desp; 
    ret &= (1 << bits) - 1;
    return ret;
}
typedef struct ENCODER {
    const char* filename;
    FILE *file;
    int dtr;
    uint32_t curr_cache, next_cache;
    uint16_t bracket_seq;
} ENCODER;

typedef struct DECODER {
    const char* filename;
    FILE *file;
    int dtr;
    uint32_t curr_cache, next_cache;
    bool end_of_file;
} DECODER;

uint16_t pullfrom_DECODER(DECODER *decoder, int bits) {
    
    if (bits > 16) bits = 16;
    uint16_t ret = 0;
    uint32_t bruh = 0;
    int bits_curr, bits_next;
    int slip = bits + decoder -> dtr -32;
    size_t els;
    if (slip > 0) {            
        if (decoder -> end_of_file) return 0; // if you go over the end of the file, return 0 always
        ret |=  (decoder -> curr_cache << slip) | (decoder -> next_cache >> (64 - slip));        
        decoder -> curr_cache = decoder -> next_cache; decoder -> dtr -= 32;
        els = fread(&decoder -> next_cache, 4, 1, decoder -> file);
        if (els == 0) decoder -> end_of_file = true;   
    } else {
        ret = get_bits_at(decoder -> curr_cache, decoder -> dtr, bits); 
    }
    printf("curr_cache: "); print_bits32(decoder -> curr_cache); printf("\n");
    printf("next_cache: "); print_bits32(decoder -> next_cache); printf("\n");
    printf("ret: "); print_bits(ret); printf("\n");
    decoder -> dtr += bits;
    return ret;
}
DECODER *new_DECODER(const char* filename) {// fread(&num2, 4, 1, file);

    /*
    FILE *file = fopen("data.bin", "rb");
    uint32_t num;
    fread(&num, 4, 1, file);
    printf("num: %u, ", num); print_bits32(num); printf("\n");
    fclose(file);
    */

    size_t els1, els2;
    DECODER *decoder = (DECODER*)malloc(sizeof(DECODER));
    decoder -> filename = filename;
    decoder -> file = fopen(filename, "rb");
    decoder -> dtr = 0; 
    decoder -> curr_cache = 0; decoder -> next_cache = 0;
    els1 = fread(&(decoder -> curr_cache), 4, 1, decoder -> file);  
    els2 = fread(&decoder -> next_cache, 4, 1, decoder -> file); 
    printf("/////////////////// CREATING DECODER..... //////////////////////////////\n");
    printf("curr_cache: "); print_bits32(decoder -> curr_cache); printf("\n");
    printf("next_cache: "); print_bits32(decoder -> next_cache); printf("\n");
    printf("////////////////////////////////////////////////////////////////////////\n");
    decoder -> end_of_file = (els2 == 0) ? true : false;
    if (els1 == 0) {
        free(decoder);
        return NULL;
    }
}
void free_DECODER(DECODER *decoder) {
    fclose(decoder -> file);
    free(decoder);
}

void pushto_ENCODER(ENCODER *encoder, uint16_t CODE, bool min_two) {
    encode_16_32(CODE, &(encoder -> curr_cache), &(encoder -> next_cache), &(encoder -> dtr), min_two);
    printf("///////////////////////////// ENCODING PRE /////////////////////////////////////\n");
    printf("curr_cache: "); print_bits32(encoder -> curr_cache); printf("\n");
    printf("next_cache: "); print_bits32(encoder -> next_cache); printf("\n");
    printf("////////////////////////////////////////////////////////////////////////////\n");

    if (encoder -> dtr >= 32) {
        printf("here\n");
        fwrite(&(encoder -> curr_cache), 4, 1, encoder -> file);
        encoder -> curr_cache = encoder -> next_cache; encoder -> next_cache = 0;
        encoder -> dtr -= 32;
    }
    printf("///////////////////////////// ENCODING AFTER /////////////////////////////////////\n");
    printf("curr_cache: "); print_bits32(encoder -> curr_cache); printf("\n");
    printf("next_cache: "); print_bits32(encoder -> next_cache); printf("\n");
    printf("////////////////////////////////////////////////////////////////////////////\n");
}
ENCODER *new_ENCODER(const char* filename) {
    ENCODER *encoder = (ENCODER*)malloc(sizeof(ENCODER));
    encoder -> filename = filename;
    encoder -> file = fopen(filename, "ab");
    encoder -> dtr = 0; encoder -> curr_cache = 0; encoder -> next_cache = 0;
    encoder -> bracket_seq = 65535;
    pushto_ENCODER(encoder, encoder -> bracket_seq, false);
    return encoder;
    // file in append mode
}
void free_ENCODER(ENCODER *encoder) {
    printf("///////////////////////////// CLOSING /////////////////////////////////////\n");
    printf("curr_cache: "); print_bits32(encoder -> curr_cache); printf("\n");
    printf("next_cache: "); print_bits32(encoder -> next_cache); printf("\n");
    printf("////////////////////////////////////////////////////////////////////////////\n");
    pushto_ENCODER(encoder, encoder -> bracket_seq, false);
    if (encoder -> dtr != 0) {
        fwrite(&(encoder -> curr_cache), 4, 1, encoder -> file);
    } 

    fclose(encoder -> file);
    free(encoder);
}


int ENCODE_DATA(const char* filename, DATA_NODE *NODES, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER) {
    ENCODER *encoder = new_ENCODER("data.bin");
    DATA_NODE *node = NODES;
    uint16_t code;
    bool encode_VAL = false;
    if (!NODES) return -1;
    while (node -> next) {
        if (encode_VAL) {
            pushto_ENCODER(encoder, node -> VAL, false);
        } else {
            if (!get_code(&code, node -> rrrrssss, CODES, VALUES, CODES_NUMBER)) return -1;  
            pushto_ENCODER(encoder, code, true);
            node = node -> next;
        } encode_VAL = !encode_VAL;
    } free_ENCODER(encoder); return 0;
}


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
    
    ENCODER *encoder = new_ENCODER("data.bin");


    pushto_ENCODER(encoder, 0, false); 


    free_ENCODER(encoder);
    /*
    FILE * file = fopen("data.bin", "wb");
    uint32_t num = 4294934527;
    fwrite(&num, 4, 1, file);
    fclose(file);*/
    /*
    file = fopen("data.bin", "rb");
    uint32_t num2;
    fread(&num2, 4, 1, file);
    printf("%u\n", num2);filename
    fclose(file);

    print_bits32(num2);

    */
    uint16_t code; DECODER *decoder = new_DECODER("data.bin");

    code = pullfrom_DECODER(decoder, 16); 
    printf("Code: %d\n", code);

    code = pullfrom_DECODER(decoder, 1); 
    printf("Code: %d\n", code);

    code = pullfrom_DECODER(decoder, 15); 
    printf("Code: %d\n", code);
    
    free_DECODER(decoder);

    


    return 0;
}







/*
    FILE *file;
    FILE *file_out;

    
    
    size_t els;

    // write
    if (true){
    file = fopen("data.bin", "wb");

    uint32_t n = 257;


    
    if (!file) return 1;
    els = fwrite(&n, sizeof(uint32_t), 1, file);
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