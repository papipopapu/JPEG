#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "core/image_compression.h"

void print_block(uint8_t *block, int size) {
    int width = 8;
    for (int i = 0; i < size; i++) {
        printf("%3u ", block[i]);
        if ((i + 1) % width == 0) printf("\n");
    }
}

/* Test basic 8x8 encode/decode roundtrip */
int test_basic_encode_decode(void) {
    printf("=== Test: Basic 8x8 Encode/Decode ===\n");
    
    uint16_t w, h;
    uint8_t *rin, *gin, *bin;
    uint8_t *rout = malloc(sizeof(uint8_t) * 64);
    uint8_t *gout = malloc(sizeof(uint8_t) * 64);
    uint8_t *bout = malloc(sizeof(uint8_t) * 64);
    
    if (!rout || !gout || !bout) {
        printf("FAIL: Memory allocation failed\n");
        return 1;
    }

    /* Initialize with random values */
    srand(42);  /* Fixed seed for reproducibility */
    for (int i = 0; i < 64; i++) {
        rout[i] = rand() % 256;
        gout[i] = rand() % 256;
        bout[i] = rand() % 256;
    }

    /* Encode */
    int ret = encode_image("test_img.bin", rout, gout, bout, 8, 8);
    if (ret != 0) {
        printf("FAIL: encode_image returned %d\n", ret);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Decode */
    ret = decode_image("test_img.bin", &rin, &gin, &bin, &w, &h);
    if (ret != 0) {
        printf("FAIL: decode_image returned %d\n", ret);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Verify dimensions */
    if (w != 8 || h != 8) {
        printf("FAIL: Dimensions mismatch. Expected 8x8, got %dx%d\n", w, h);
        free(rin); free(gin); free(bin);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Cleanup */
    free(rin); free(gin); free(bin);
    free(rout); free(gout); free(bout);
    
    /* Remove test file */
    remove("test_img.bin");
    
    printf("PASS\n\n");
    return 0;
}

/* Test larger image (multiple blocks) */
int test_multi_block_image(void) {
    printf("=== Test: Multi-Block 16x16 Image ===\n");
    
    uint16_t width = 16, height = 16;
    size_t size = width * height;
    uint16_t w_out, h_out;
    uint8_t *rin, *gin, *bin;
    uint8_t *rout = malloc(sizeof(uint8_t) * size);
    uint8_t *gout = malloc(sizeof(uint8_t) * size);
    uint8_t *bout = malloc(sizeof(uint8_t) * size);
    
    if (!rout || !gout || !bout) {
        printf("FAIL: Memory allocation failed\n");
        return 1;
    }

    /* Initialize with gradient pattern */
    for (size_t i = 0; i < size; i++) {
        rout[i] = (uint8_t)(i % 256);
        gout[i] = (uint8_t)((i * 2) % 256);
        bout[i] = (uint8_t)((i * 3) % 256);
    }

    /* Encode */
    int ret = encode_image("test_multi.bin", rout, gout, bout, width, height);
    if (ret != 0) {
        printf("FAIL: encode_image returned %d\n", ret);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Decode */
    ret = decode_image("test_multi.bin", &rin, &gin, &bin, &w_out, &h_out);
    if (ret != 0) {
        printf("FAIL: decode_image returned %d\n", ret);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Verify dimensions */
    if (w_out != width || h_out != height) {
        printf("FAIL: Dimensions mismatch. Expected %dx%d, got %dx%d\n", 
               width, height, w_out, h_out);
        free(rin); free(gin); free(bin);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Cleanup */
    free(rin); free(gin); free(bin);
    free(rout); free(gout); free(bout);
    remove("test_multi.bin");
    
    printf("PASS\n\n");
    return 0;
}

/* Test solid color image (uniform values) */
int test_solid_color(void) {
    printf("=== Test: Solid Color Image ===\n");
    
    uint16_t width = 8, height = 8;
    size_t size = width * height;
    uint16_t w_out, h_out;
    uint8_t *rin, *gin, *bin;
    uint8_t *rout = malloc(sizeof(uint8_t) * size);
    uint8_t *gout = malloc(sizeof(uint8_t) * size);
    uint8_t *bout = malloc(sizeof(uint8_t) * size);
    
    if (!rout || !gout || !bout) {
        printf("FAIL: Memory allocation failed\n");
        return 1;
    }

    /* Initialize with solid color (128, 128, 128) */
    for (size_t i = 0; i < size; i++) {
        rout[i] = 128;
        gout[i] = 128;
        bout[i] = 128;
    }

    /* Encode */
    int ret = encode_image("test_solid.bin", rout, gout, bout, width, height);
    if (ret != 0) {
        printf("FAIL: encode_image returned %d\n", ret);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Decode */
    ret = decode_image("test_solid.bin", &rin, &gin, &bin, &w_out, &h_out);
    if (ret != 0) {
        printf("FAIL: decode_image returned %d\n", ret);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* Verify dimensions */
    if (w_out != width || h_out != height) {
        printf("FAIL: Dimensions mismatch. Expected %dx%d, got %dx%d\n", 
               width, height, w_out, h_out);
        free(rin); free(gin); free(bin);
        free(rout); free(gout); free(bout);
        return 1;
    }

    /* For solid color, expect very close values after roundtrip */
    int max_diff = 0;
    for (size_t i = 0; i < size; i++) {
        int diff_r = abs((int)rin[i] - 128);
        int diff_g = abs((int)gin[i] - 128);
        int diff_b = abs((int)bin[i] - 128);
        if (diff_r > max_diff) max_diff = diff_r;
        if (diff_g > max_diff) max_diff = diff_g;
        if (diff_b > max_diff) max_diff = diff_b;
    }
    
    /* Solid color should have minimal distortion (tolerance for DCT/quantization) */
    if (max_diff > 20) {
        printf("FAIL: Solid color distortion too high: max_diff=%d\n", max_diff);
        free(rin); free(gin); free(bin);
        free(rout); free(gout); free(bout);
        return 1;
    }
    printf("  Max difference from original: %d (acceptable)\n", max_diff);

    /* Cleanup */
    free(rin); free(gin); free(bin);
    free(rout); free(gout); free(bout);
    remove("test_solid.bin");
    
    printf("PASS\n\n");
    return 0;
}

/* Test error handling for invalid filename */
int test_decode_invalid_file(void) {
    printf("=== Test: Decode Invalid File ===\n");
    
    uint16_t w, h;
    uint8_t *r, *g, *b;
    
    /* Try to decode non-existent file */
    int ret = decode_image("nonexistent_file.bin", &r, &g, &b, &w, &h);
    if (ret == 0) {
        printf("FAIL: decode_image should fail for non-existent file\n");
        free(r); free(g); free(b);
        return 1;
    }
    
    /* Try to decode with NULL filename */
    ret = decode_image(NULL, &r, &g, &b, &w, &h);
    if (ret == 0) {
        printf("FAIL: decode_image should fail for NULL filename\n");
        free(r); free(g); free(b);
        return 1;
    }
    
    /* Try to decode with empty filename */
    ret = decode_image("", &r, &g, &b, &w, &h);
    if (ret == 0) {
        printf("FAIL: decode_image should fail for empty filename\n");
        free(r); free(g); free(b);
        return 1;
    }
    
    printf("PASS\n\n");
    return 0;
}

int main(void) {
    int failed = 0;
    
    printf("JPEG Compression Library Tests\n");
    printf("==============================\n\n");
    
    failed += test_basic_encode_decode();
    failed += test_multi_block_image();
    failed += test_solid_color();
    failed += test_decode_invalid_file();
    
    printf("==============================\n");
    if (failed == 0) {
        printf("All tests PASSED!\n");
        return 0;
    } else {
        printf("%d test(s) FAILED!\n", failed);
        return 1;
    }
}