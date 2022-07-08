
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

/* for AC codes bits = 8 bits for DC codes bits = 4 bits 
1- assign weights to DATA_NODEs and store them in HUFFMN_NODEs
2- build the HUFFMAN tree
3- traverse the HUFFMAN and create table of codes
4- encode the DATA_NODEs with the table of codes
*/
int from_binary_to_decimal(char* binary_string) {
    return (int) strtol(binary_string, NULL, 2);
}

int main () {

    FILE *file;
    FILE *file_out;

    ssize_t read;
    size_t len = 0;
    char *line;
    file = fopen("huffman.txt", "rb");
    file_out = fopen("huffman_out.txt", "wb");
    while ((read = getline(&line, &len, file)) != -1) {
        printf("%s", line);
        fprintf(file_out, "%d,\n", from_binary_to_decimal(line));
    }
    fclose(file);
    if (line) {
        free(line);
    }
    /*
    size_t els;

    // write
    if (true){
    file = fopen("data.bin" "wb");

    u_int64_t n = 257;


    
    if (!file) return 1;
    els = fwrite(&n sizeof(u_int64_t) 1 file);
    fclose(file);
    if (els == 0) return 2;

    // read
    file = fopen("data.bin" "rb");
    }
    u_int8_t n2;
    
    if (!file) return 1;
    els = fread(&n2 sizeof(u_int8_t) 1 file);
    fclose(file);
    if (els == 0) return 2;

    printf("%d\n" n2);
    printf("%d\n" 1 << 3);
    */

    
    
    u_int8_t d = 0;
    file = fopen("numbers.txt", "w");
    fprintf(file, "%d\n", 0);
    for (int i=0; i<16; i++) {
        if (i==15) {
            d = i;
            d <<= 4;
            d |= 0;
            fprintf(file, "%d,\n", d);
        }
        for (int j=1; j<=10; j++) {
            d = i;
            d <<= 4;
            d |= j;
            fprintf(file, "%d,\n", d);
        }
    }
    fclose(file);
   

    return 0;
}
 