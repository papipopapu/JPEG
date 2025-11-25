#include "image_compression.h"


void get_block(uint8_t *slice, uint8_t *UINT8_BLOCK, uint16_t IMG_WIDTH, uint16_t IMG_HEIGHT, int I0, int J0) {
    /*
    Extract a block of size 8 from IMAGE at position (I0, J0). If the block overflows the image, it is filled
    with an approximation based on near values to reach 8 * 8 pixels.
        * IMAGE: the image to extract the block from.
        * UINT8_BLOCK: the block to fill outside the image.
        * IMG_WIDTH: the width of the image.
        * IMG_HEIGHT: the height of the image.
        * I0: the x coordinate of upper left corner of the block.
        * J0: the y coordinate of upper left corner of the block.
    */  

    int i, j, disti, distj;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            
            if (j + J0 >= IMG_WIDTH && i + I0 >= IMG_HEIGHT) {
                distj = 2+ J0 + j - IMG_WIDTH; 
                disti = 2+ I0 + i - IMG_HEIGHT;
                UINT8_BLOCK[i * 8 + j] = (UINT8_BLOCK[(i-1) * 8 + j] + UINT8_BLOCK[i * 8 + (j-1)]) / (disti + distj);
            }
            else if (i + I0 > IMG_HEIGHT) {
                disti = 2+ I0 + i - IMG_HEIGHT;
                UINT8_BLOCK[i * 8 + j] = UINT8_BLOCK[(i-1) * 8 + j] / disti;
            }
            else if (j + J0 > IMG_WIDTH) {
                distj = 2+ J0 + j - IMG_WIDTH;
                UINT8_BLOCK[i * 8 + j] = UINT8_BLOCK[i * 8 + (j-1)] / distj;
            }
            else if (j + J0 == IMG_WIDTH) {
                UINT8_BLOCK[i * 8 + j] = slice[(I0 + i) * IMG_WIDTH + (J0 + j-1)] / 2.;
            }       
            else if (i + I0 == IMG_HEIGHT) {
                UINT8_BLOCK[i * 8 + j] = slice[(I0 + i-1) * IMG_WIDTH + (J0 + j)] / 2.;
            }
            else {
                UINT8_BLOCK[i * 8 + j] = slice[(I0 + i) * IMG_WIDTH + (J0 + j)];
            }
        }
    }
}
void put_block(uint8_t *slice, uint8_t *UINT8_BLOCK, uint16_t IMG_WIDTH, uint16_t IMG_HEIGHT, int I0, int J0) {
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            if (j+J0 < IMG_WIDTH && i+I0 < IMG_HEIGHT) {
                slice[(I0 + i)*IMG_WIDTH + J0 + j] = UINT8_BLOCK[i * 8 + j];
            }
        }
    }
}

void block_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK) {
    /*
    Obtains the discrete cosine transform of the given BLOCK of pixeks, into the FLOAT_BLOCK, both of size 8 * 8.
    Args:
        * UINT8_BLOCK: input block
        * FLOAT_BLOCK: output block, icontains dct transformation          
    */
    int i, j, k, l;
    float ai, aj, temp, cte = 2./8;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            ai = i == 0 ? M_SQRT1_2 : 1;
            aj = j == 0 ? M_SQRT1_2 : 1;  
            temp = 0;      
            for (k = 0; k < 8; k++) {
                for (l = 0; l < 8; l++) {
                    temp +=  (UINT8_BLOCK[k * 8 + l] - 128) * cos(M_PI*i*0.5*(2.*k+1.)/8) * cos(M_PI*j*0.5*(2.*l+1.)/8);
                }
            } 
            FLOAT_BLOCK[i * 8 + j] = cte * ai * aj * temp; 
        }
    }
}

void block_inv_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK) 
{
    /*
    Obtains the inverse discrete cosine transform of the given BLOCK of pixels, into the FLOAT_BLOCK, both of size 8 * 8.
    Args:
        * UINT8_BLOCK: output block
        * FLOAT_BLOCK: input block, contains the dct transform           
    */
    int u, v, x, y;
    float au, av, temp, cte = 2./8;
    for (x = 0;  x < 8; x++) {
        for (y = 0; y < 8; y++) {
            temp = 0;      
            for (u = 0; u < 8; u++) {
                for (v = 0; v < 8; v++) {
                    au = u == 0 ? M_SQRT1_2 : 1;
                    av = v == 0 ? M_SQRT1_2 : 1;  
                    temp += au * av * FLOAT_BLOCK[u * 8 + v] * cos(M_PI*u*0.5*(2.*x+1.)/8) * cos(M_PI*v*0.5*(2.*y+1.)/8);
                }
            } 
            int val = (int)round(cte * temp) + 128;
            if (val < 0) val = 0;
            if (val > 255) val = 255;
            UINT8_BLOCK[x * 8 + y] = (uint8_t)val; 
        }
    }
}

void general_dct(uint8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_WIDTH, size_t BLOCK_HEIGHT) {
    /*
    Obtains the discrete cosine transform of the given chunk of pixeks, into the FLOAT_BLOCK.
    Used when a whole BLOCK does not fit, here only part of the BLOCK varaible's memory will be used.
    Args:
        * UINT8_BLOCK: input block. The block is assumed to be of size BLOCK_WIDTH * BLOCK_HEIGHT.
        * FLOAT_BLOCK: output block
        * BLOCK_WIDTH: width of the block
        * BLOCK_HEIGHT: height of the block
    */
    int i, j, k, l;
    float ai, aj, temp, cte = 2./sqrt(BLOCK_WIDTH * BLOCK_HEIGHT);
    for (i = 0; i < BLOCK_HEIGHT; i++) {
        for (j = 0; j < BLOCK_WIDTH; j++) {
            ai = i == 0 ? M_SQRT1_2 : 1;
            aj = j == 0 ? M_SQRT1_2 : 1;  
            temp = 0;      
            for (k = 0; k < BLOCK_HEIGHT; k++) {
                for (l = 0; l < BLOCK_WIDTH; l++) {
                    temp +=  UINT8_BLOCK[k * BLOCK_WIDTH + l] * cos(M_PI*i*0.5*(2.*k+1.)/BLOCK_HEIGHT) * cos(M_PI*j*0.5*(2.*l+1.)/BLOCK_WIDTH);
                }
            } 
            FLOAT_BLOCK[i * BLOCK_WIDTH + j] = cte * ai * aj * temp; 
        }
    }
}

void block_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK) {
    /*
    Quantizies the FLOAT_BLOCK after a discrete cosine transform using the QUANT_MAT.
    Args:
        * QUANT_MAT: quantization matrix
        * INT16_BLOCK: output block
        * FLOAT_BLOCK: input block
    */
    int i;
    for (i = 0; i < 8 * 8; i++) {  
        INT16_BLOCK[i] = round(FLOAT_BLOCK[i] / QUANT_MAT[i]);  
    }
}

void block_inv_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK) {
    /*
    Obtains the inverse quantization of the INT16_BLOCK using the QUANT_MAT.
    Args:
        * QUANT_MAT: quantization matrix
        * FLOAT_BLOCK: output block 
        * INT16_BLOCK: input block
    */
    int i;
    for (i = 0; i < 64; i++) {  
        FLOAT_BLOCK[i] = INT16_BLOCK[i] * QUANT_MAT[i];  
    }
}

void block_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX) {
    /*
    Reorders the INT8_BLOCK into INT16_SEQUENCE after a quantization.
    Args:
        * INT16_SEQUENCE: output block
        * INT16_BLOCK: input block
    */
    for (int i = 0; i < 64; i++) {
        INT16_SEQUENCE[i] = INT16_BLOCK[SERIAL_IDX[i]];
    }
}


void block_inv_serialize(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, const uint8_t *SERIAL_IDX) {
    /*
    Undoes the serialize reordering of the INT16_SEQUENCE into the INT16_BLOCK.
    Args:
        * INT16_SEQUENCE: input block
        * INT16_BLOCK: output block
    */
    for (int i = 0; i < 64; i++) {
        INT16_BLOCK[SERIAL_IDX[i]] = INT16_SEQUENCE[i];
    }
}

void DATA_PACKET_pack(DATA_PACKET *data, int16_t VAL, uint8_t zeros) {
    /*
    Pack a value into a DATA_PACKET for JPEG encoding.
    
    JPEG coefficient encoding:
    - Category (ssss) = minimum bits needed to represent |VAL|
    - For category n, values range from -2^n+1 to -2^(n-1) and 2^(n-1) to 2^n-1
    - Positive values: write the value directly (n bits)
    - Negative values: write (value - 1) in one's complement form (n bits)
    
    Examples:
    - VAL = 0:  ssss = 0, no additional bits
    - VAL = 1:  ssss = 1, bits = 1
    - VAL = -1: ssss = 1, bits = 0
    - VAL = 2:  ssss = 2, bits = 10
    - VAL = -2: ssss = 2, bits = 01
    - VAL = 3:  ssss = 2, bits = 11
    - VAL = -3: ssss = 2, bits = 00
    */
    int min_bits;
    int16_t abs_val = (VAL < 0) ? -VAL : VAL;
    
    min_bits = (abs_val == 0) ? 0 : min_bits_abs(abs_val);
    
    data->rrrrssss = zeros;
    data->rrrrssss <<= 4;
    data->rrrrssss |= (uint8_t)min_bits;
    data->VAL_bits = min_bits;
    
    if (VAL > 0) {
        data->VAL = (uint16_t)VAL;
    } else if (VAL < 0) {
        // One's complement: flip all bits of |VAL|, keeping only min_bits
        data->VAL = (uint16_t)(VAL - 1) & ((1 << min_bits) - 1);
    } else {
        data->VAL = 0;
    }
}

bool DATA_PACKET_encode(DATA_PACKET *data, const uint16_t *CODES, const uint8_t *VALUES, const uint8_t *LENGTHS, size_t N_CODES) {
    int i;
    for (i = 0; i < N_CODES; i++) {
        if (VALUES[i] == data -> rrrrssss) {
            data -> rs_code = CODES[i];
            data -> rs_code_bits = LENGTHS[i];
            return true;
        }
    }
    printf("Not found matching code for %d\n", data -> rrrrssss);
    return false;
}

int block_encode(OUTSTREAM* out, int16_t *INT16_SEQUENCE, int16_t *PREV_DC,
 const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS, 
 const uint16_t *AC_CODES, const uint8_t *AC_VALUES, const uint8_t *AC_LENGTHS) {
    int i, zeros = 0; int16_t val = INT16_SEQUENCE[0];
    int last_nonzero_idx = 0;
    
    // Find last non-zero AC coefficient
    for (i = 63; i >= 1; i--) {
        if (INT16_SEQUENCE[i] != 0) {
            last_nonzero_idx = i;
            break;
        }
    }

    DATA_PACKET data;
    DATA_PACKET_pack(&data, val - *PREV_DC, 0);
    if (!DATA_PACKET_encode(&data, DC_CODES, DC_VALUES, DC_LENGTHS, 12)) return -1;
    OUTSTREAM_push(out, data.rs_code, data.rs_code_bits);
    OUTSTREAM_push(out, data.VAL, data.VAL_bits);
    *PREV_DC = val;
    
    for (i = 1; i <= last_nonzero_idx; i++) {
        val = INT16_SEQUENCE[i];
        if (val == 0 && zeros < 15) {
            zeros++;
        } else {
            DATA_PACKET_pack(&data, val, zeros);
            if (!DATA_PACKET_encode(&data, AC_CODES, AC_VALUES, AC_LENGTHS, 162)) return -1;
            OUTSTREAM_push(out, data.rs_code, data.rs_code_bits);
            OUTSTREAM_push(out, data.VAL, data.VAL_bits);
            zeros = 0;
        }
    }
    
    // Write EOB if there are trailing zeros (or if all AC coefficients are zero)
    if (last_nonzero_idx < 63) {
        OUTSTREAM_push(out, AC_CODES[0], AC_LENGTHS[0]);
    }
    return 0;
}

int block_decode(INSTREAM* in, int16_t *INT16_SEQUENCE, int16_t* PREV_DC,
 const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS,
 const uint16_t *AC_CODES, const uint8_t *AC_VALUES,const uint8_t *AC_LENGTHS) {

    uint16_t val;
    uint8_t rrrrssss;
    int rrrr, ssss, idx=1;
    bool eob = false;
    
    // Initialize sequence to zeros
    for (int i = 0; i < 64; i++) INT16_SEQUENCE[i] = 0;

    // Decode DC coefficient
    if(!search_codes(in, &rrrrssss, DC_CODES, DC_VALUES, DC_LENGTHS, 12)) {printf("DC Code not found!\n"); return -1;}
    decode_data(in, rrrrssss, &ssss, &rrrr, &val);
    if (rrrr != 0) return -1; // integrity check: DC should have run=0
    
    // Decode the DC differential value
    int16_t dc_diff;
    if (ssss == 0) {
        dc_diff = 0;
    } else if (val & (1 << (ssss - 1))) {
        dc_diff = (int16_t)val;
    } else {
        dc_diff = (int16_t)val - (1 << ssss) + 1;
    }
    INT16_SEQUENCE[0] = dc_diff + *PREV_DC;
    *PREV_DC = INT16_SEQUENCE[0];

    // Decode AC coefficients
    while (!eob && idx < 64) {
        if(!search_codes(in, &rrrrssss, AC_CODES, AC_VALUES, AC_LENGTHS, 162)){printf("AC Code not found at idx %d!\n", idx); return -1;}
        decode_data(in, rrrrssss, &ssss, &rrrr, &val);
        eob = write_data(INT16_SEQUENCE, false, idx, ssss, rrrr, val);
        idx++; idx += rrrr;
    }    
    return 0;
}


int min_bits_abs(int16_t n) {  
    if (n == 0) return 1;
    n = (n<0) ? -n : n;
    int i;
    int count = 0;
    for (i = 15; i >= 0; i--) {
        if ((n >> i) & 1) return 16 - count;
        count++;
    }
    return 0;
}

int min_bits_code(uint16_t n) {  
    if (n == 0 || n == 1) return 2;
    int i, count = 0;
    for (i = 15; i >= 0; i--) {
        if ((n >> i) & 1) return 16 - count;
        count++;
    }
    return 0;
}



bool search_codes(INSTREAM *in, uint8_t *rrrrssss, const uint16_t *CODES, const uint8_t *VALUES, const uint8_t *LENGTHS, size_t CODES_NUMBER) {
    /* Checks if there are any matches of any amount of crecent digits of the compare base inside the CODES provided. */
    uint16_t code = 0, pull = 0;
    int bits = 2, i;
    INSTREAM_pull(in, &code, 2);
    
    while(bits <= 16) {
        for (i = 0; i < CODES_NUMBER; i++) {
            if ((code == CODES[i]) && (bits == LENGTHS[i])) {
                *rrrrssss = VALUES[i];
                return true;
            }
        }
        if (bits == 16) break;  // Already checked 16-bit codes, don't read more
        INSTREAM_pull(in, &pull, 1);
        code = (code << 1) | pull;
        bits++;
    }
    return false;
}

void decode_data(INSTREAM *in, uint8_t rrrrssss, int *ssss, int *rrrr, uint16_t *val) {
    /* Decodes the data packet. */
    *rrrr = rrrrssss >> 4;
    *ssss = rrrrssss & 0x0F;
    
    *val = 0; INSTREAM_pull(in, val, *ssss);
    // printf("Pulling %d bits TO GET VAL: ", *ssss); print_16bits(*val); printf("\n");
    

}

bool write_data(int16_t *INT16_SEQUENCE, bool is_dc, int idx, int ssss, int rrrr, uint16_t val) {
    /*
    Decode and write a value to the sequence.
    
    JPEG coefficient decoding:
    - If MSB of val (in ssss bits) is 1, value is positive: true_val = val
    - If MSB of val is 0, value is negative: true_val = val - 2^ssss + 1
    
    Examples:
    - ssss = 1, val = 1 -> true_val = 1
    - ssss = 1, val = 0 -> true_val = 0 - 2 + 1 = -1
    - ssss = 2, val = 10 (2) -> true_val = 2
    - ssss = 2, val = 01 (1) -> true_val = 1 - 4 + 1 = -2
    - ssss = 2, val = 11 (3) -> true_val = 3
    - ssss = 2, val = 00 (0) -> true_val = 0 - 4 + 1 = -3
    */
    
    // Check for EOB (End of Block)
    if (!is_dc && rrrr == 0 && ssss == 0) {
        // EOB - remaining coefficients are zero (already initialized)
        return true;
    }

    // Write run of zeros before the value
    for (int k = 0; k < rrrr && (idx + k) < 64; k++) {
        INT16_SEQUENCE[idx + k] = 0;
    }
    idx += rrrr;
    
    if (idx >= 64) return false;

    // Decode and write the value
    if (ssss == 0) {
        INT16_SEQUENCE[idx] = 0;
        return false;
    }
    
    // Check MSB to determine sign
    int16_t true_val;
    if (val & (1 << (ssss - 1))) {
        // MSB is 1: positive value
        true_val = (int16_t)val;
    } else {
        // MSB is 0: negative value (one's complement)
        true_val = (int16_t)val - (1 << ssss) + 1;
    }
    
    INT16_SEQUENCE[idx] = true_val;
    return false;
}
