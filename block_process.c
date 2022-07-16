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
    int i, j, disti, distj;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            if (j+J0 < IMG_WIDTH && i+I0 < IMG_HEIGHT) {
                slice[(I0 + i)*IMG_WIDTH + J0 + j] = UINT8_BLOCK[i * 8 + j];
            }
        }
    }
}
void block_downsample420(uint8_t *UINT8_BLOCK)
{
    /*
    Downsample a block of size 8 with a ratio of 4:2:0.
    Args:
        * UINT8_BLOCK: the block to downsample.
    */
    int i;
    for (i = 0; i < 8; i++) {
        UINT8_BLOCK[i] = 2 * round(UINT8_BLOCK[i] / 2.);  
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
    Obtains the inversse discrete cosine transform of the given BLOCK of pixeks, into the FLOAT_BLOCK, both of size 8 * 8.
    Args:
        * UINT8_BLOCK: output block
        * FLOAT_BLOCK: inùt block, contains the dct transform           
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
            UINT8_BLOCK[x * 8 + y] = round(cte * temp) + 128; 
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
    for (int i = 0; i < 8 * 8; i++) {
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
    for (int i = 0; i < 8 * 8; i++) {
        INT16_SEQUENCE[i] = INT16_BLOCK[SERIAL_IDX[8 * 8 - i - 1]];
    }
}

void DATA_PACKET_pack(DATA_PACKET *data, int16_t VAL, uint8_t zeros) {
    uint8_t is_neg = 0; 
    int min_bits; 
    if (VAL < 0) {VAL = -VAL;  is_neg = 1;}
    min_bits = (VAL == 0) ? 0 : min_bits_abs(VAL);
    data -> rrrrssss = zeros; 
    data -> rrrrssss <<= 4; 
    data -> rrrrssss |= min_bits; 
    data -> VAL = VAL;
    data -> VAL &= (1 << (min_bits-1))-1; 
    data -> VAL |= is_neg << (min_bits-1); 
    data -> VAL_bits = min_bits;

    //printf("Packed rrrrssss: %d, bits: ", data -> rrrrssss); print_ubits(data -> rrrrssss); printf("\n");
    
    // 0 size -> read 0 (val ~+-0      ) -> { read 0-1  (+) append sign bit } absurd, we know its 0 
    // 1 size -> read 1 (val ~+-1      ) -> { read 1-1  (+) append sign bit }
    // ...
    // n size -> read n (val ~+-2^(n-1)) -> { read n-1  (+) append sign bit }

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
void print_matrix(int16_t *seq) {
    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 8; j++) {
            printf("%d ", seq[i*8+j]);
        }
        printf("\n");
    }
}
int block_encode(OUTSTREAM* out, int16_t *INT16_SEQUENCE, int16_t *PREV_DC,
 const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS, 
 const uint16_t *AC_CODES, const uint8_t *AC_VALUES, const uint8_t *AC_LENGTHS) {
    int i, zeros = 0; int16_t val = INT16_SEQUENCE[0];
    printf("Recieved sequence: \n"); print_matrix(INT16_SEQUENCE); printf("\n");

    DATA_PACKET data;
    DATA_PACKET_pack(&data, val - *PREV_DC, 0);
    if (!DATA_PACKET_encode(&data, DC_CODES, DC_VALUES, DC_LENGTHS, 12)) return -1;
    OUTSTREAM_push(out, data.rs_code, data.rs_code_bits);
    OUTSTREAM_push(out, data.VAL, data.VAL_bits);
    *PREV_DC = val;
    
    printf("//////////////////////////////////////////////////////////////////\n");
    printf("rrrrssss, bits="); print_16bits(data.rrrrssss); printf("\n");
    printf("rs code, length=%d, bits=", data.rs_code_bits); print_16bits(data.rs_code); printf("\n");
    printf("value=%d, length=%d, bits=", val, data.VAL_bits); print_16bits(data.VAL); printf("\n");
    

    
    for (i = 1; i < 64; i++) {
        val = INT16_SEQUENCE[i];
        if (val == 0 && zeros < 15) {
            zeros++;
        } else {
            DATA_PACKET_pack(&data, val, zeros);
            if (!DATA_PACKET_encode(&data, AC_CODES, AC_VALUES, AC_LENGTHS, 162)) return -1;
            OUTSTREAM_push(out, data.rs_code, data.rs_code_bits);
            OUTSTREAM_push(out, data.VAL, data.VAL_bits);
            zeros = 0;
    
    printf("//////////////////////////////////////////////////////////////////\n");
    printf("rrrrssss, bits="); print_16bits(data.rrrrssss); printf("\n");
    printf("rs code, length=%d, bits=", data.rs_code_bits); print_ubits(data.rs_code); printf("\n");
    printf("value=%d, length=%d, bits=", val, data.VAL_bits); print_16bits(data.VAL); printf("\n");
    
        }
    }
    
    printf("//////////////////////////////////////////////////////////////////e\n");
    printf("eob\n");

    OUTSTREAM_push(out, AC_CODES[0], AC_LENGTHS[0]);
    return 0;
}

int block_decode(INSTREAM* in, int16_t *INT16_SEQUENCE, int16_t* PREV_DC,
 const uint16_t *DC_CODES, const uint8_t *DC_VALUES, const uint8_t *DC_LENGTHS,
 const uint16_t *AC_CODES, const uint8_t *AC_VALUES,const uint8_t *AC_LENGTHS) {

    uint16_t val;
    uint8_t rrrrssss;
    int rrrr, ssss, idx=1;
    bool eob = false;
    

    if(!search_codes(in, &rrrrssss, DC_CODES, DC_VALUES, DC_LENGTHS, 12)) {printf("Code not found!\n"); return -1;}
    decode_data(in, rrrrssss, &ssss, &rrrr, &val);
    if (rrrr != 0) return -1; // some integrity checking, non dc read (run != 0)
    write_data(INT16_SEQUENCE, true, 0, ssss, 0, val+*PREV_DC);
    *PREV_DC = INT16_SEQUENCE[0];


    while (!eob && idx < 64) {
        if(!search_codes(in, &rrrrssss, AC_CODES, AC_VALUES, AC_LENGTHS, 162)){printf("Code not found!\n"); return -1;};
        decode_data(in, rrrrssss, &ssss, &rrrr, &val);
        eob = write_data(INT16_SEQUENCE, false, idx, ssss, rrrr, val);
        idx++; idx+= rrrr;
    }
    printf("Idx: %d, Eob: %d\n", idx, eob);
    
    return eob ? 0 : 1;
}
//////////////////////////



void print_f(float *matrix) {
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            printf("%f ", matrix[i * 8 + j]);
        }
        printf("\n");
    }
}
void print_ubits(uint8_t n)
{
    int i;
    for (i = 7; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
}
void print_bits(int8_t n)
{
    int i;
    for (i = 7; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
}
void print_16bits(int16_t n)
{
    int i;
    for (i = 15; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
}
void print_32bits(uint32_t n)
{
    int i;
    for (i = 31; i >= 0; i--) {
        printf("%d", (n >> i) & 1);
    }
}



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
    int bits = 2, i, nigga;
    INSTREAM_pull(in, &code, 2);
    printf("LOoking for next rrrrsss\n");
    while(bits < 16) {
        printf("Checking bits=%d, code: ", bits); print_16bits(code); printf("\n");
        for (i = 0; i < CODES_NUMBER; i++) {
            //printf("Comparing with: bits=%d, code=", (int)LENGTHS[i]); print_16bits(CODES[i]); printf("\n");
            if ((code == CODES[i]) && (bits == LENGTHS[i])) {
                *rrrrssss = VALUES[i];
               
                printf("Found code! i=%d, rrrrssss=", i); print_ubits(*rrrrssss); printf("\n");
                return true;
            }
        }
        INSTREAM_pull(in, &pull, 1);
        code = (code << 1) | pull;
        bits++;
    }
    printf("Not found!\n");
    return false;
}

void decode_data(INSTREAM *in, uint8_t rrrrssss, int *ssss, int *rrrr, uint16_t *val) {
    /* Decodes the data packet. */
    *rrrr = rrrrssss >> 4;
    *ssss = rrrrssss & 0b1111;
    
    *val = 0; INSTREAM_pull(in, val, *ssss);
    printf("Pulling %d bits TO GET VAL: ", *ssss); print_16bits(*val); printf("\n");
    

}

bool write_data(int16_t *INT16_SEQUENCE, bool is_dc, int idx, int ssss, int rrrr, uint16_t val) {
    // check if eob
    //printf("Run: %d\n", rrrr);
    if (!is_dc && rrrr == 0 && ssss == 0) {
        while(idx < 64) {INT16_SEQUENCE[idx++] = 0;}
        return true;
    }

    // write zeros
    int k = 0, safety_counter = 0;
    while(k<rrrr) {INT16_SEQUENCE[idx+k] = 0; k++; safety_counter++; if (safety_counter > 15) {
        printf("Exceeded safety"); return false;}}
    idx += rrrr;

    // then, we write the value

    // special cases
    if (ssss == 0) {INT16_SEQUENCE[idx] = 0;  return false;}
    if (ssss == 1) {INT16_SEQUENCE[idx] = (val == 1) ? -1 : 1; return false;}

    // normal cases
    int i; bool is_neg; int16_t true_val;
    is_neg = (val & (1 << (ssss - 1))) != 0; // check sign bit
    val |= (1 << (ssss - 1)); // add the missing 1
    if (is_neg) {true_val = -val;} else {true_val = val;}
    INT16_SEQUENCE[idx] = true_val;
    printf("true val pulled: %d\n", true_val);
    return false;
 
}
