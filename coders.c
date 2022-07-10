#include "image_compression.h"
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

void encode_to_cache(uint16_t CODE, uint32_t *curr_cache, uint32_t *next_cache, int *dtr, bool min_two) {
    //*next_cache = [00000000...];  curr_cache = [{prev_codes}0000000...]
    uint32_t bigger = CODE;
    int bits = min_bits(CODE); if (min_two) bits += (CODE == 1 | CODE == 0) ? 1 : 0;
    if (*dtr + bits > 32) {
        *curr_cache |= bigger >> bits - 32 + *dtr; 
        *next_cache |= bigger << 64 - bits - *dtr;
    } else if (*dtr + bits == 32) {
        *curr_cache |= bigger;
    } else {
        *curr_cache |= bigger << 32 - bits - *dtr;
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
    size_t els;
    if (decoder -> dtr + bits > 32) {            
        if (decoder -> end_of_file) return 0;
        bits_curr = 32 - decoder -> dtr; bits_next = bits - bits_curr;
        ret |=  (get_bits_at(decoder -> curr_cache, decoder -> dtr     , bits_curr) << bits_next)
            |   (get_bits_at(decoder -> next_cache, decoder -> dtr - 32, bits_next)             );
        decoder -> curr_cache = decoder -> next_cache;
        els = fread(&decoder -> next_cache, 4, 1, decoder -> file);
        decoder -> dtr -= 32;
        if (els == 0) {decoder -> end_of_file = true; printf("Reached end of file\n");}      
    } else {
        ret = get_bits_at(decoder -> curr_cache, decoder -> dtr, bits); 
    }
    decoder -> dtr += bits;
    
    return ret;
}
DECODER *new_DECODER(const char* filename) {
    size_t els1, els2;
    DECODER *decoder = (DECODER*)malloc(sizeof(DECODER));
    decoder -> filename = filename;
    decoder -> file = fopen(filename, "rb");
    decoder -> dtr = 0; 
    els1 = fread(&decoder -> curr_cache, 4, 1, decoder -> file);  
    els2 = fread(&decoder -> next_cache, 4, 1, decoder -> file); 
    decoder -> end_of_file = (els2 == 0) ? true : false;
    return (els1 == 0) ? NULL : decoder;
}
void free_DECODER(DECODER *decoder) {
    fclose(decoder -> file);
    free(decoder);
}

void pushto_ENCODER(ENCODER *encoder, uint16_t CODE, bool min_two) {
    encode_to_cache(CODE, &(encoder -> curr_cache), &(encoder -> next_cache), &(encoder -> dtr), min_two);
    if (encoder -> dtr >= 32) {
        fwrite(&(encoder -> curr_cache), 4, 1, encoder -> file);
        encoder -> curr_cache = encoder -> next_cache; encoder -> next_cache = 0;
        encoder -> dtr -= 32;
    }
}
ENCODER *new_ENCODER(const char* filename) {
    ENCODER *encoder = (ENCODER*)malloc(sizeof(ENCODER));
    encoder -> filename = filename;
    encoder -> file = fopen(filename, "wb");
    encoder -> dtr = 0; encoder -> curr_cache = 0; encoder -> next_cache = 0;
    encoder -> bracket_seq = 65535;
    pushto_ENCODER(encoder, encoder -> bracket_seq, false);
    return encoder;
    // file in append mode
}
void free_ENCODER(ENCODER *encoder) {
    pushto_ENCODER(encoder, encoder -> bracket_seq, false);
    if (encoder -> curr_cache != 0) {
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

int DECODE_DATA(const char* filename, DATA_NODE **NODES, const uint16_t *CODES, const uint8_t *VALUES, size_t CODES_NUMBER) {
    DECODER *decoder = new_DECODER("data.bin");
    DATA_NODE *node = NULL;
    uint16_t code;
    bool decode_VAL = false;
    if (!NODES) return -1;
    while (true) {
        if (decode_VAL) {
            node = new_DATA_NODE(node, pullfrom_DECODER(decoder, 16), NULL);
        } else {
            if (!get_code(&code, pullfrom_DECODER(decoder, 16), CODES, VALUES, CODES_NUMBER)) return -1;
            node = new_DATA_NODE(node, code, NULL);
            decode_VAL = !decode_VAL;
        } if (decoder -> end_of_file) break;
    } free_DECODER(decoder);
    *NODES = node; return 0;
}