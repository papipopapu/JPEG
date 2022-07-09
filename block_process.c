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
uint8_t min_bits_abs(int16_t n) {
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

uint8_t min_bits(int16_t n) {
    // min bits to hold the value of int16 //
    int i;
    uint8_t count = 0;
    for (i = 15; i >= 0; i--) {
        if ((n >> i) & 1) return 16 - count;
        count++;
    }
    return 0;
}


void block_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, bool IS_FIRST) {

    /*
    Packs the data from one blcck into the AC-DC linked list.
    Args:
        * INT16_SEQUENCE: sequence of ints  of 1 block.
        * AC_DATA_NODES: pointer new  pointer to the latest node of a linked list of AC PACKED_NODEs.
        * DC_DATA_NODES: pointer new  pointer to the latest node of a linked list of DC PACKED_NODEs.
        * IS_FIRST: true if this is the first block of the image.
    */
    int i;
    uint8_t zeros = 0;
    int16_t val, prev_DC_val = 0, temp; 

    DATA_NODE *prev_DC, *curr_DC, *prev_AC, *curr_AC;

    if (IS_FIRST) { // SPECIAL CASE
        prev_AC = NULL; prev_DC = NULL;
    } else {
        prev_AC = *AC_DATA_NODES; prev_DC = *DC_DATA_NODES;
    }

    curr_DC = new_DATA_NODE(); pack_DATA_NODE(curr_DC, 0, INT16_SEQUENCE[0]); // DC value
    connect_DATA_NODE(&prev_DC, &curr_DC, DC_DATA_NODES);

    for (i = 1; i < 8 * 8; i++) { 
        val = INT16_SEQUENCE[i];
        if (val == 0 && zeros < 15) { //  running
            zeros++; 
        } else {
            curr_AC = new_DATA_NODE(); pack_DATA_NODE(curr_AC, zeros, val); // ACs and 16 zeros runs
            connect_DATA_NODE(&prev_AC, &curr_AC, AC_DATA_NODES);
            zeros = 0;
        }   
    }
    
    curr_AC = new_DATA_NODE(); pack_DATA_NODE(curr_AC, 0, 0); // EOB
    prev_AC -> next = curr_AC;    
}
void blocks_pack(int16_t *INT16_SEQUENCE, DATA_NODE **AC_DATA_NODES, DATA_NODE **DC_DATA_NODES, size_t BLOCK_NUMBER) {
    /*
    Packs the data from the swhole image sequence into a linked list of structs. 
    Args:
        * INT16_SEQUENCE: sequence of ints for whole image (until EOB).
        * AC_DATA_NODES: pointer new  pointer to the head of a linked list of AC PACKED_NODEs.
        * DC_DATA_NODES: pointer new  pointer to the head of a linked list of DC PACKED_NODEs.
        * BLOCK_NUMBER: number of blocks in the image.
    */
    int i;
    if (BLOCK_NUMBER == 0) return;
    block_pack(INT16_SEQUENCE, AC_DATA_NODES, DC_DATA_NODES, true);
    for (i = 1; i < BLOCK_NUMBER; i ++) { 
          block_pack(INT16_SEQUENCE + i * 64, AC_DATA_NODES, DC_DATA_NODES, false);
    }

    // size = 0 -> val = 0
    // (15, 0) (0) is the usual compressed token
    // (0, 0) (0) EOB collides with the unlikely scenario of an ending 0 preceded by a filled AC, but when decompressing, we gon be
    // keeping track of the dimensions of the original image, so we can know when to stop (i.e. stop = all data read && EOB reached, else throw ERROR)
}

//////////////////////////
void print_ui8(uint8_t *matrix) {
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            printf("%d ", matrix[i * 8 + j]);
        }
        printf("\n");
    }
}
void print_i8(int8_t *matrix) {
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            printf("%d ", matrix[i * 8 + j]);
        }
        printf("\n");
    }
}
void print_i16(int16_t *matrix) {
    int i, j;
    for (i = 0; i < 8; i++) {
        for (j = 0; j < 8; j++) {
            printf("%d ", matrix[i * 8 + j]);
        }
        printf("\n");
    }
}
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
void print_list(DATA_NODE* head)
{
    DATA_NODE* curr = head;
    printf("\n");
    while (curr != NULL) {
        printf("zeros_bitsVAL: ");
        print_ubits(curr->zeros_bitsVAL);
        printf(", VAL bits: ");
        print_16bits(curr->VAL);
        printf("\n");
        
        curr = curr->next;
    }
    printf("\n");
}
/////////////////////////////
void block_process_one(bool isY, uint8_t *UINT8_BLOCK, DATA_NODE **AC_HEAD, DATA_NODE **DC_HEAD)
{
    size_t TOTAL_BITSIZE;
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nInitial block:\n");
    print_ui8(UINT8_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Downsample
    if (!isY) {
        block_downsample420(UINT8_BLOCK);
    } 
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nDownsampled:\n");
    print_ui8(UINT8_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // DCT
    float *FLOAT_BLOCK = malloc(sizeof(float) * 64);
    block_dct(UINT8_BLOCK, FLOAT_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nAfter dct:\n");
    print_f(FLOAT_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // inverse DCT
    block_inv_dct(UINT8_BLOCK, FLOAT_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nReconstructed after idct:\n");
    print_ui8(UINT8_BLOCK);
    block_dct(UINT8_BLOCK, FLOAT_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Quantize
    int16_t *INT16_BLOCK = malloc(sizeof(int16_t) * 64);
    if (isY) {
        block_quantize(LUMINANCE_QUANT, INT16_BLOCK, FLOAT_BLOCK);
    } else {
        block_quantize(CHROMINANCE_QUANT, INT16_BLOCK, FLOAT_BLOCK);
    }
    block_quantize(CHROMINANCE_QUANT, INT16_BLOCK, FLOAT_BLOCK);
    free(FLOAT_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nQuantized:\n");
    print_i16(INT16_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    // serialize
    int16_t *INT16_SEQUENCE = malloc(sizeof(int16_t) * 64);
    block_serialize(INT16_BLOCK, INT16_SEQUENCE, ZIGZAG_IDX);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nSequenced:\n");
    print_i16(INT16_SEQUENCE);
    ////////////////////////////////////////////////////////////////////////////////////////////////
     // Inverse serialize
    block_inv_serialize(INT16_BLOCK, INT16_SEQUENCE, ZIGZAG_IDX);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\n De-sequenced:\n");
    print_i16(INT16_BLOCK);
    block_serialize(INT16_BLOCK, INT16_SEQUENCE, ZIGZAG_IDX);
    free(INT16_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Pack
    blocks_pack(INT16_SEQUENCE, AC_HEAD, DC_HEAD, 1);
    free(INT16_SEQUENCE);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nPacked DCs:\n");
    print_list(*DC_HEAD);
    printf("\nPacked ACs:\n");
    print_list(*AC_HEAD);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    //printf("\nPrevious size: 512, new size: %ld\n", TOTAL_BITSIZE);
}
/*
uint16 bits_to_write
uint16 extra_bits
uint8 rem_bits_num
while (node->next)
 bits_to_write = 
 extra_bits = (node -> next) ? next -> getcode(node -> next -> val) : 
 */
// {uint8_block} -> dct -> {float_block} -> quantize -> {int8_block} -> serialize_reorder -> {int8_sequence} -> huffman -> {huffman}
void blocks_encode(FILE* file, DATA_NODE *AC_DATA_NODES, DATA_NODE *DC_DATA_NODES,
    const uint16_t *DC_NEWCODES, const uint16_t *DC_OLDCODES, const uint16_t *AC_NEWCODES, const uint16_t *AC_OLDCODES,
    size_t BLOCK_NUMBER, size_t CODES_NUMBER) {
    /* Encodes a set of blocks into a file.
    Args:
     * file: The file to write to.
     * AC_DATA_NODES: The AC data nodes.
     * DC_DATA_NODES: The DC data nodes.
     * DC_NEWCODES: The new DC codes.
     * DC_OLDCODES: The old DC codes.
     * AC_NEWCODES: The new AC codes.
     * AC_OLDCODES: The old AC codes.
     * 8: The block dimension.
     * BLOCK_NUMBER: The number of blocks.
     * CODES_NUMBER: The number of codes.
     */
    // 16 ob + DC + 16ob + AC + 16ob + fill bits
    size_t encoded_bits = 0, i;
    uint16_t ob = 65535;
    fopen(file, "a");
    fwrite(ob, sizeof(uint16_t), 1, file);
    


    

    
}

