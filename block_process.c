#include "image_compression.h"
static size_t bhs;
static size_t ahs;

void get_block(uint8_t *IMAGE, uint8_t *UINT8_BLOCK, size_t IMG_WIDTH, size_t IMG_HEIGHT, size_t I0, size_t J0) {
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
                UINT8_BLOCK[i * 8 + j] = IMAGE[(I0 + i) * IMG_WIDTH + (J0 + j-1)] / 2.;
            }       
            else if (i + I0 == IMG_HEIGHT) {
                UINT8_BLOCK[i * 8 + j] = IMAGE[(I0 + i-1) * IMG_WIDTH + (J0 + j)] / 2.;
            }
            else {
                UINT8_BLOCK[i * 8 + j] = IMAGE[(I0 + i) * IMG_WIDTH + (J0 + j)];
            }
        }
    }
}
void block_rgb_to_yCbCr(uint8_t *r_to_y, uint8_t *g_to_Cb, uint8_t *b_to_Cr)
{
    /*
    Translates rgb values to yCbCr values.
        * r_to_y: the block containing r values, and output for y.
        * g_to_Cb: the block containing g values, and output for Cb.
        * b_to_Cr: the block containing b values, and output for Cr.
    */
   int i;
   uint8_t R, G, B;
   for (i = 0; i < 8 * 8; i++) {
        R = r_to_y[i];
        G = g_to_Cb[i];
        B = b_to_Cr[i];
        r_to_y[i] = (R * 0.299 + G * 0.587 + B * 0.114);
        g_to_Cb[i] = (R * -0.168736 + G * -0.331264 + B * 0.5 + 128);
        b_to_Cr[i] = (R * 0.5 + G * -0.418688 + B * -0.081312 + 128);
   }
}
void block_yCbCt_to_rgb(uint8_t *y_to_r, uint8_t *Cb_to_g, uint8_t *Cr_to_b)
{
    /*
    Translates yCbCr values to rgb values.
        * y_to_r: the block containing y values, and output for r.
        * Cb_to_g: the block containing Cb values, and output for g.
        * Cr_to_b: the block containing Cr values, and output for b.
    */
   int i;
   uint8_t Y, Cb, Cr;
   for (i = 0; i < 8 * 8; i++) {
        Y = y_to_r[i];
        Cb = Cb_to_g[i];
        Cr = Cr_to_b[i];
        y_to_r[i] = (Y + 1.402 * (Cr - 128));
        Cb_to_g[i] = (Y - 0.344136 * (Cb - 128) - 0.714136 * (Cr - 128));
        Cr_to_b[i] = (Y + 1.772 * (Cb - 128));
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
    
    // 0 size -> read 0 (val ~+-0      ) -> { read 0-1  (+) append sign bit } absurd, we know its 0 
    // 1 size -> read 1 (val ~+-1      ) -> { read 1-1  (+) append sign bit }
    // ...
    // n size -> read n (val ~+-2^(n-1)) -> { read n-1  (+) append sign bit }

}

bool DATA_PACKET_encode(DATA_PACKET *data, const uint16_t *CODES, const uint8_t *VALUES, size_t N_CODES) {
    int i;
    for (i = 0; i < N_CODES; i++) {
        if (VALUES[i] == data -> rrrrssss) {
            data -> rs_code = CODES[i];
            data -> rs_code_bits = min_bits_code(data -> rs_code);
            return true;
        }
    }
    return false;
}

int block_encode(OUTSTREAM* out, int16_t *INT16_SEQUENCE, int16_t PREV_DC,
 const uint16_t *DC_CODES, const uint8_t *DC_VALUES,
 const uint16_t *AC_CODES, const uint8_t *AC_VALUES) {
    int i, zeros = 0; int16_t val = INT16_SEQUENCE[0];
    DATA_PACKET data;
    DATA_PACKET_pack(&data, val - PREV_DC, 0);
    if (!DATA_PACKET_encode(&data, DC_CODES, DC_VALUES, 12)) return -1;
    OUTSTREAM_push(out, data.rs_code, data.rs_code_bits);
    OUTSTREAM_push(out, data.VAL, data.VAL_bits);

    printf("//////////////////////////////////////////////////////////////////\n");
    printf("rrrrssss, bits="); print_16bits(data.rrrrssss); printf("\n");
    printf("rs code, length=%d, bits=", data.rs_code_bits); print_16bits(data.rs_code); printf("\n");
    printf("value=%d, length=%d, bits=", val, data.VAL_bits); print_16bits(data.VAL); printf("\n");


    
    for (i = 1; i < 64; i++) {
        val = INT16_SEQUENCE[i];
        if (val == 0 && zeros < 16) {
            zeros++;
        } else {
            DATA_PACKET_pack(&data, val, zeros);
            if (!DATA_PACKET_encode(&data, AC_CODES, AC_VALUES, 162)) return -1;
            OUTSTREAM_push(out, data.rs_code, data.rs_code_bits);
            OUTSTREAM_push(out, data.VAL, data.VAL_bits);
            zeros = 0;
    printf("//////////////////////////////////////////////////////////////////\n");
    printf("rrrrssss, bits="); print_16bits(data.rrrrssss); printf("\n");
    printf("rs code, length=%d, bits=", data.rs_code_bits); print_16bits(data.rs_code); printf("\n");
    printf("value=%d, length=%d, bits=", val, data.VAL_bits); print_16bits(data.VAL); printf("\n");

        }
    }
    DATA_PACKET_pack(&data, val, zeros);
    if (!DATA_PACKET_encode(&data, AC_CODES, AC_VALUES, 162)) return -1;
    OUTSTREAM_push(out, data.rs_code, data.rs_code_bits);
    OUTSTREAM_push(out, data.VAL, data.VAL_bits);

    printf("//////////////////////////////////////////////////////////////////\n");
    printf("rrrrssss, bits="); print_16bits(data.rrrrssss); printf("\n");
    printf("rs code, length=%d, bits=", data.rs_code_bits); print_16bits(data.rs_code); printf("\n");
    printf("value=%d, length=%d, bits=", val, data.VAL_bits); print_16bits(data.VAL); printf("\n");

    OUTSTREAM_push(out, 0, 8);
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
    if (n == 0 | n == 1) return 2;
    int i, count = 0;
    for (i = 15; i >= 0; i--) {
        if ((n >> i) & 1) return 16 - count;
        count++;
    }
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
