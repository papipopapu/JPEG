#include "image_compression.h"
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
bool get_code(uint16_t *code, uint8_t rrrrssss, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER) {
    int i;
    for (i=0; i < CODES_NUMBER; i++) {
        if (rrrrssss == VALUES[i]) {*code = CODES[i]; return true;}
    } return false;
}
void encode_to_cache(uint16_t CODE, uint32_t *curr_cache, uint32_t *next_cache, int *dtr, bool HUFFMAN_MODE) {
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
    decoder -> dtr += bits;
    return ret;
}
DECODER *new_DECODER(const char* filename) {// fread(&num2, 4, 1, file);
    size_t els1, els2;
    DECODER *decoder = (DECODER*)malloc(sizeof(DECODER));
    decoder -> filename = filename;
    decoder -> file = fopen(filename, "rb");
    decoder -> dtr = 0; 
    els1 = fread(&(decoder -> curr_cache), 4, 1, decoder -> file);  
    els2 = fread(&decoder -> next_cache, 4, 1, decoder -> file); 
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
    encode_to_cache(CODE, &(encoder -> curr_cache), &(encoder -> next_cache), &(encoder -> dtr), min_two);
    if (encoder -> dtr >= 32) {
        printf("here\n");
        fwrite(&(encoder -> curr_cache), 4, 1, encoder -> file);
        encoder -> curr_cache = encoder -> next_cache; encoder -> next_cache = 0;
        encoder -> dtr -= 32;
    }
}
ENCODER *new_ENCODER(const char* filename) {
    ENCODER *encoder = (ENCODER*)malloc(sizeof(ENCODER));
    encoder -> filename = filename;
    encoder -> file = fopen(filename, "ab");
    encoder -> dtr = 0; encoder -> curr_cache = 0; encoder -> next_cache = 0;
    encoder -> bracket_seq = 65535;
    pushto_ENCODER(encoder, encoder -> bracket_seq, false);
    return encoder;
}
void free_ENCODER(ENCODER *encoder) {
    pushto_ENCODER(encoder, encoder -> bracket_seq, false);
    if (encoder -> dtr != 0) {
        fwrite(&(encoder -> curr_cache), 4, 1, encoder -> file);
    } 
    fclose(encoder -> file);
    free(encoder);
}


int ENCODE_DATA(const char* filename, DATA_NODE *NODES, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER) {
    ENCODER *encoder = new_ENCODER(filename);
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