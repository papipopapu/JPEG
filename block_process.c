#include "image_compression.h"


void get_block(u_int8_t *IMAGE, u_int8_t *UINT8_BLOCK, size_t BLOCK_SIZE, size_t IMG_WIDTH, size_t IMG_HEIGHT, size_t I0, size_t J0) {
    /*
    Extract a block of size BLOCK_SIZE from IMAGE at position (I0, J0). If the block overflows the image, it is filled
    with an approximation based on near values to reach BLOCK_SIZE * BLOCK_SIZE pixels.
        * IMAGE: the image to extract the block from.
        * UINT8_BLOCK: the block to fill outside the image.
        * BLOCK_SIZE: the size of the block.
        * IMG_WIDTH: the width of the image.
        * IMG_HEIGHT: the height of the image.
        * I0: the x coordinate of upper left corner of the block.
        * J0: the y coordinate of upper left corner of the block.
    */  

    int i, j, disti, distj;
    for (i = 0; i < BLOCK_SIZE; i++) {
        for (j = 0; j < BLOCK_SIZE; j++) {
            
            if (j + J0 >= IMG_WIDTH && i + I0 >= IMG_HEIGHT) {
                distj = 2+ J0 + j - IMG_WIDTH; 
                disti = 2+ I0 + i - IMG_HEIGHT;
                UINT8_BLOCK[i * BLOCK_SIZE + j] = (UINT8_BLOCK[(i-1) * BLOCK_SIZE + j] + UINT8_BLOCK[i * BLOCK_SIZE + (j-1)]) / (disti + distj);
            }
            else if (i + I0 > IMG_HEIGHT) {
                disti = 2+ I0 + i - IMG_HEIGHT;
                UINT8_BLOCK[i * BLOCK_SIZE + j] = UINT8_BLOCK[(i-1) * BLOCK_SIZE + j] / disti;
            }
            else if (j + J0 > IMG_WIDTH) {
                distj = 2+ J0 + j - IMG_WIDTH;
                UINT8_BLOCK[i * BLOCK_SIZE + j] = UINT8_BLOCK[i * BLOCK_SIZE + (j-1)] / distj;
            }
            else if (j + J0 == IMG_WIDTH) {
                UINT8_BLOCK[i * BLOCK_SIZE + j] = IMAGE[(I0 + i) * IMG_WIDTH + (J0 + j-1)] / 2.;
            }       
            else if (i + I0 == IMG_HEIGHT) {
                UINT8_BLOCK[i * BLOCK_SIZE + j] = IMAGE[(I0 + i-1) * IMG_WIDTH + (J0 + j)] / 2.;
            }
            else {
                UINT8_BLOCK[i * BLOCK_SIZE + j] = IMAGE[(I0 + i) * IMG_WIDTH + (J0 + j)];
            }
        }
    }
}
void bloc_rgb_to_yCbCr(u_int8_t *r_to_y, u_int8_t *g_to_Cb, u_int8_t *b_to_Cr, size_t BLOCK_DIM)
{
    /*
    Translates rgb values to yCbCr values.
        * r_to_y: the block containing r values, and output for y.
        * g_to_Cb: the block containing g values, and output for Cb.
        * b_to_Cr: the block containing b values, and output for Cr.
        * BLOCK_DIM: the size of the block.
    */
   int i;
   u_int8_t R, G, B;
   for (i = 0; i < BLOCK_DIM * BLOCK_DIM; i++) {
        R = r_to_y[i];
        G = g_to_Cb[i];
        B = b_to_Cr[i];
        r_to_y[i] = (R * 0.299 + G * 0.587 + B * 0.114);
        g_to_Cb[i] = (R * -0.168736 + G * -0.331264 + B * 0.5 + 128);
        b_to_Cr[i] = (R * 0.5 + G * -0.418688 + B * -0.081312 + 128);
   }
}
void block_downsample420(u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM)
{
    /*
    Downsample a block of size BLOCK_SIZE with a ratio of 4:2:0.
    Args:
        * UINT8_BLOCK: the block to downsample.
        * BLOCK_DIM: block dimension, 8 for 8x8 blocks, 16 for 16x16 blocks, etc.   
    */
    int i;
    for (i = 0; i < BLOCK_DIM; i++) {
        UINT8_BLOCK[i] = 2 * round(UINT8_BLOCK[i] / 2.);  
    }
}
void block_dct(u_int8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_DIM) {
    /*
    Obtains the discrete cosine transform of the given BLOCK of pixeks, into the FLOAT_BLOCK, both of size BLOCK_DIM * BLOCK_DIM.
    Args:
        * UINT8_BLOCK: input block, its a pointer to an isolated chunk of memory outside the original image.
        * FLOAT_BLOCK: output block, its a pointer to an isolated chunk of memory outside the original image. 
        * BLOCK_DIM: block dimension, 8 for 8x8 blocks, 16 for 16x16 blocks, etc.     
          
    */
    int i, j, k, l;
    float ai, aj, temp, cte = 2./BLOCK_DIM;
    for (i = 0; i < BLOCK_DIM; i++) {
        for (j = 0; j < BLOCK_DIM; j++) {
            ai = i == 0 ? M_SQRT1_2 : 1;
            aj = j == 0 ? M_SQRT1_2 : 1;  
            temp = 0;      
            for (k = 0; k < BLOCK_DIM; k++) {
                for (l = 0; l < BLOCK_DIM; l++) {
                    temp +=  (UINT8_BLOCK[k * BLOCK_DIM + l] - 128) * cos(M_PI*i*0.5*(2.*k+1.)/BLOCK_DIM) * cos(M_PI*j*0.5*(2.*l+1.)/BLOCK_DIM);
                }
            } 
            FLOAT_BLOCK[i * BLOCK_DIM + j] = cte * ai * aj * temp; 
        }
    }
}

void block_idct(float *FLOAT_BLOCK, u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM)
{
    /*
    Obtains the inversse discrete cosine transform of the given BLOCK of pixeks, into the FLOAT_BLOCK, both of size BLOCK_DIM * BLOCK_DIM.
    Args:
        * UINT8_BLOCK: output block, its a pointer to an isolated chunk of memory outside the original image.
        * FLOAT_BLOCK: inùt block, contains the dct transform 
        * BLOCK_DIM: block dimension, 8 for 8x8 blocks, 16 for 16x16 blocks, etc.     
          
    */
    int u, v, x, y;
    float au, av, temp, cte = 2./BLOCK_DIM;
    for (x = 0;  x < BLOCK_DIM; x++) {
        for (y = 0; y < BLOCK_DIM; y++) {
            temp = 0;      
            for (u = 0; u < BLOCK_DIM; u++) {
                for (v = 0; v < BLOCK_DIM; v++) {
                    au = u == 0 ? M_SQRT1_2 : 1;
                    av = v == 0 ? M_SQRT1_2 : 1;  
                    temp += au * av * FLOAT_BLOCK[u * BLOCK_DIM + v] * cos(M_PI*u*0.5*(2.*x+1.)/BLOCK_DIM) * cos(M_PI*v*0.5*(2.*y+1.)/BLOCK_DIM);
                }
            } 
            UINT8_BLOCK[x * BLOCK_DIM + y] = round(cte * temp) + 128; 
        }
    }
}

void general_dct(u_int8_t *UINT8_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_WIDTH, size_t BLOCK_HEIGHT) {
    /*
    Obtains the discrete cosine transform of the given chunk of pixeks, into the FLOAT_BLOCK.
    Used when a whole BLOCK does not fit, here only part of the BLOCK varaible's memory will be used.
    Args:
        * UINT8_BLOCK: input block, its a pointer to an isolated chunk of memory outside the original image. 
                           The block is assumed to be of size BLOCK_WIDTH * BLOCK_HEIGHT.
        * FLOAT_BLOCK: output block, its a pointer to an isolated chunk of memory outside the original image.
        * BLOCK_DIM: block dimension, 8 for 8x8 blocks, 16 for 16x16 blocks, etc.     
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

void block_quantize(const float *QUANT_MAT, int16_t *INT16_BLOCK, float *FLOAT_BLOCK, size_t BLOCK_DIM) {
    /*
    Quantizies the FLOAT_BLOCK after a discrete cosine transform using the QUANT_MAT.
    Args:
        * QUANT_MAT: quantization matrix
        * INT16_BLOCK: output block, its a pointer to an isolated chunk of memory outside image of int8 type
        * FLOAT_BLOCK: input block, its a pointer to an isolated chunk of memory outside image of float type
        * BLOCK_DIM: block dimension, 8 for 8x8 blocks, 16 for 16x16 blocks, etc.
    */
    int i;
    for (i = 0; i < BLOCK_DIM * BLOCK_DIM; i++) {  
        INT16_BLOCK[i] = round(FLOAT_BLOCK[i] / QUANT_MAT[i]);  
    }
}

void block_zigzag(int16_t *INT16_BLOCK, int16_t *INT16_SEQUENCE, size_t BLOCK_DIM) {
    /*
    Reorders the INT8_BLOCK  into INT16_SEQUENCE after a quantization.
    Args:
        * INT16_SEQUENCE: output block, its a pointer to an isolated chunk of memory outside image 
        * INT16_BLOCK: output block, its a pointer to an isolated chunk of memory outside image 
        * BLOCK_DIM: block dimension, 8 for 8x8 blocks, 16 for 16x16 blocks, etc.
    */
    int i = 0, j = 0, p = 0;
    INT16_SEQUENCE[p++] = INT16_BLOCK[i * BLOCK_DIM + j];
    while (p < BLOCK_DIM * BLOCK_DIM)
	{
		if (j < BLOCK_DIM - 1)
		{
			j++;
		}
		else
		{
			i++;
		}

		while (i < BLOCK_DIM && j >= 0)
		{
			INT16_SEQUENCE[p++] = INT16_BLOCK[i * BLOCK_DIM + j];
			i++;
			j--;
		}
		i--;
		j++;

		if (i < BLOCK_DIM - 1)
		{
			i++;
		}
		else
		{
			j++;
		}

		while (i >= 0 && j < BLOCK_DIM)
		{
			INT16_SEQUENCE[p++] = INT16_BLOCK[i * BLOCK_DIM + j];
			i--;
			j++;
		}
		i++;
		j--;
	}
}

u_int8_t min_bits_abs(int16_t n) {
    // min bits to hold the abolute value of int16 //
    int i;
    if (n < 0) n = -n;
    u_int8_t count = 0;
    for (i = 15; i >= 0; i--) {
        if ((n >> i) & 1) return 16 - count;
        count++;
    }
    return 0;
}


void block_pack(int16_t *INT16_SEQUENCE, DATA_NODE *PACKED_BLOCK_HEAD, size_t BLOCK_DIM) {
    /*
    Packs the data from the sequence into a linked list of structs. 
    Args:
        * INT16_SEQUENCE: sequence of ints.
        * PACKED_BLOCK_HEAD: pointer to the head of a linked list of PACKED_NODEs.
        * BLOCK_DIM: block dimension, 8 for 8x8 blocks, 16 for 16x16 blocks, etc.
    */
    int i;
    u_int8_t zeros = 0;
    int16_t val; 
     // DC component
    DATA_NODE *prev = PACKED_BLOCK_HEAD, *next = NULL;
    pack_DATA_NODE(prev, 0, INT16_SEQUENCE[0]);

    for (i = 1; i < BLOCK_DIM * BLOCK_DIM - 1; i++) { // AC components
        val = INT16_SEQUENCE[i];
        if (val == 0 && zeros < 15) {
            zeros++; 
        } else {
            next = new_DATA_NODE();
            pack_DATA_NODE(next, zeros, val);
            prev -> next = next;
            prev = next;
            zeros = 0;
        }   
    }
    next = new_DATA_NODE(); // AC at the end
    pack_DATA_NODE(next, zeros, 0);
    prev -> next = next;
}

//////////////////////////
void print_ui8(u_int8_t *matrix, int n, int m) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%d ", matrix[i * m + j]);
        }
        printf("\n");
    }
}
void print_i8(int8_t *matrix, int n, int m) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%d ", matrix[i * m + j]);
        }
        printf("\n");
    }
}
void print_i16(int16_t *matrix, int n, int m) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%d ", matrix[i * m + j]);
        }
        printf("\n");
    }
}
void print_f(float *matrix, int n, int m) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            printf("%f ", matrix[i * m + j]);
        }
        printf("\n");
    }
}
void print_ubits(u_int8_t n)
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
        printf("VAL: %d, zeros_bitsVAL: ", curr->VAL);
        print_ubits(curr->zeros_bitsVAL);
        printf(", VAL bits: ");
        print_16bits(curr->VAL);
        printf("\n");
        
        curr = curr->next;
    }
    printf("\n");
}
/////////////////////////////
void block_process_one(bool isY, u_int8_t *UINT8_BLOCK, size_t BLOCK_DIM, DATA_NODE *PACKED_BLOCK_HEAD)
{
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nInitial block:\n");
    print_ui8(UINT8_BLOCK, BLOCK_DIM, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Downsample
    if (!isY) {
        block_downsample420(UINT8_BLOCK, BLOCK_DIM);
    } 
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nDownsampled:\n");
    print_ui8(UINT8_BLOCK, BLOCK_DIM, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // DCT
    float *FLOAT_BLOCK = malloc(sizeof(float) * BLOCK_DIM * BLOCK_DIM);
    block_dct(UINT8_BLOCK, FLOAT_BLOCK, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nAfter dct:\n");
    print_f(FLOAT_BLOCK, BLOCK_DIM, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // inverse DCT
    block_idct(FLOAT_BLOCK, UINT8_BLOCK, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nReconstructed after idct:\n");
    print_ui8(UINT8_BLOCK, BLOCK_DIM, BLOCK_DIM);
    block_dct(UINT8_BLOCK, FLOAT_BLOCK, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Quantize
    int16_t *INT16_BLOCK = malloc(sizeof(int16_t) * BLOCK_DIM * BLOCK_DIM);
    if (isY) {
        block_quantize(LUMINANCE_QUANT_MATRIX_8_8, INT16_BLOCK, FLOAT_BLOCK, BLOCK_DIM);
    } else {
        block_quantize(CHROMINANCE_QUANT_MATRIX_8_8, INT16_BLOCK, FLOAT_BLOCK, BLOCK_DIM);
    }
    block_quantize(LUMINANCE_QUANT_MATRIX_8_8, INT16_BLOCK, FLOAT_BLOCK, BLOCK_DIM);
    free(FLOAT_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nQuantized:\n");
    print_i16(INT16_BLOCK, BLOCK_DIM, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Zigzag
    int16_t *INT16_SEQUENCE = malloc(sizeof(int16_t) * BLOCK_DIM * BLOCK_DIM);
    block_zigzag(INT16_BLOCK, INT16_SEQUENCE, BLOCK_DIM);
    free(INT16_BLOCK);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nSequenced:\n");
    print_i16(INT16_SEQUENCE, BLOCK_DIM, BLOCK_DIM);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Pack
    block_pack(INT16_SEQUENCE, PACKED_BLOCK_HEAD, BLOCK_DIM);
    free(INT16_SEQUENCE);
    ////////////////////////////////////////////////////////////////////////////////////////////////
    printf("\nPacked list:\n");
    print_list(PACKED_BLOCK_HEAD);
    ////////////////////////////////////////////////////////////////////////////////////////////////
}

// {uint8_block} -> dct -> {float_block} -> quantize -> {int8_block} -> zigzag_reorder -> {int8_sequence} -> huffman -> {huffman}

