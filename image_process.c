#include "image_compression.h"

RGB_IMAGE *new_RGB_IMAGE(uint16_t width, uint16_t height) {
    RGB_IMAGE *img = malloc(sizeof(RGB_IMAGE));
    img->WIDTH = width;
    img->HEIGHT = height;
    img->r = malloc(sizeof(uint8_t) * width * height);
    img->g = malloc(sizeof(uint8_t) * width * height);
    img->b = malloc(sizeof(uint8_t) * width * height);
    return img;
}
void delete_RGB_IMAGE(RGB_IMAGE *img) {
    free(img->r);
    free(img->g);
    free(img->b);
    free(img);
}
void downsample_420(uint8_t * Cb, uint8_t *Cr, size_t N) { // must be ycbcr
    for (int i = 0; i < N; i++) {
        Cb[i] = 2*round(Cb[i]/2.0);
        Cr[i] = 2*round(Cr[i]/2.0);
    }
}

void image_rgb_to_yCbCr(uint8_t *r_to_y, uint8_t *g_to_Cb, uint8_t *b_to_Cr)
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
void image_yCbCr_to_rgb(uint8_t *y_to_r, uint8_t *Cb_to_g, uint8_t *Cr_to_b)
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
int encode_slice(OUTSTREAM *out, uint8_t *slice, uint16_t WIDTH, uint16_t HEIGHT, bool is_luminance) {
    int I, J;
    uint8_t *block = malloc(sizeof(uint8_t) * 64);
    int16_t *sequence = malloc(sizeof(int16_t) * 64), *int16_block = malloc(sizeof(int16_t) * 64);
    float * float_block = malloc(sizeof(float) * 64);

    int16_t PREV_DC = 0;
    for (I=0; I<HEIGHT; I+=8) {
        for (J=0; J<WIDTH; J+=8) {
            get_block(slice, block, WIDTH, HEIGHT, I, J);
            block_dct(block, float_block);
            if (is_luminance) {
                block_quantize(LUMINANCE_QUANT , int16_block, float_block );
                block_serialize(int16_block, sequence,  ZIGZAG_IDX);
                if (block_encode(out, sequence, &PREV_DC, DC_LUMINANCE_CODES, DC_VALUES, DC_LUMINANCE_LENGTHS,
                                                      AC_LUMINANCE_CODES, AC_VALUES, AC_LUMINANCE_LENGTHS) != 0) {printf("Error encoding..."); return -1;}
            } else {
                block_quantize(CHROMINANCE_QUANT , int16_block, float_block );
                block_serialize(int16_block, sequence,  ZIGZAG_IDX);
                if (block_encode(out, sequence, &PREV_DC, DC_CHROMINANCE_CODES, DC_VALUES, DC_CHROMINANCE_LENGTHS,
                                                      AC_CHROMINANCE_CODES, AC_VALUES, AC_CHROMINANCE_LENGTHS) != 0) {printf("Error encoding..."); return -1;};
            }  

        }
    }
    free(block); free(sequence); free(int16_block); free(float_block);
    return 0;
}
int encode_image(char *filename, uint8_t *r, uint8_t *g, uint8_t *b,  uint16_t width, uint16_t height) {
    OUTSTREAM *out = new_OUTSTREAM(filename, 8);
    OUTSTREAM_push(out, width, 16); OUTSTREAM_push(out, height, 16);
    image_rgb_to_yCbCr(r, g, b);
    downsample_420(g, b, width * height);
    // r = y, g = Cb, b = Cr
    // downsample
    encode_slice(out, r, width, height, true);
    encode_slice(out, g, width, height, false);
    encode_slice(out, b, width, height, false);
    delete_OUTSTREAM(out);
    return 0;
}
int decode_slice(INSTREAM *in, uint8_t *slice, uint16_t WIDTH, uint16_t HEIGHT, bool is_luminance) {
    int I, J;
    uint8_t *block = malloc(sizeof(uint8_t) * 64);
    int16_t *sequence = malloc(sizeof(int16_t) * 64), *int16_block = malloc(sizeof(int16_t) * 64);
    float * float_block = malloc(sizeof(float) * 64);
    int16_t PREV_DC = 0;
    for (I=0; I<HEIGHT; I+=8) {
        for (J=0; J<WIDTH; J+=8) {
            if (is_luminance) {
            if (block_decode(in, sequence, &PREV_DC, DC_LUMINANCE_CODES, DC_VALUES, DC_LUMINANCE_LENGTHS,
                                                     AC_LUMINANCE_CODES, AC_VALUES, AC_LUMINANCE_LENGTHS) != 0) {printf("Error decoding..."); return -1;}
                
                block_inv_serialize(int16_block, sequence, ZIGZAG_IDX);
                block_inv_quantize(LUMINANCE_QUANT, int16_block, float_block);

            } else {
            if (block_decode(in, sequence, &PREV_DC, DC_CHROMINANCE_CODES, DC_VALUES, DC_CHROMINANCE_LENGTHS,
                                                     AC_CHROMINANCE_CODES, AC_VALUES, AC_CHROMINANCE_LENGTHS) != 0) {printf("Error decoding..."); return -1;};
                block_inv_serialize(int16_block, sequence, ZIGZAG_IDX);
                block_inv_quantize(CHROMINANCE_QUANT, int16_block, float_block);
            }
            block_inv_dct(block, float_block);
            put_block(slice, block, WIDTH, HEIGHT, I, J);
        }
    }
    free(block); free(sequence); free(int16_block); free(float_block);
    return 0;
}
int decode_image(char *filename, uint8_t **r, uint8_t **g, uint8_t **b, uint16_t *width, uint16_t *height) {
    if (filename == NULL || strlen(filename) == 0){
        printf("Error in arguments provided\n");
        return -1;
    }
    
    INSTREAM *in = new_INSTREAM(filename, 8);
    INSTREAM_pull(in, width, 16); INSTREAM_pull(in, height, 16);

    *r = malloc(sizeof(uint8_t) * (*width) * (*height));
    *g = malloc(sizeof(uint8_t) * (*width) * (*height));
    *b = malloc(sizeof(uint8_t) * (*width) * (*height));

    decode_slice(in, *r, *width, *height, true);
    decode_slice(in, *g, *width, *height, false);
    decode_slice(in, *b, *width, *height, false);
    image_yCbCr_to_rgb(*r, *g, *b);
    delete_INSTREAM(in);
    return 0;
}

