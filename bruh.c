#include "stdio.h"
#include "stdlib.h"
#include <string.h>

int main() {
    int length[162];
    // count the number of characters on each line of a file, and store them in length[]
    FILE *file = fopen("numbers.txt", "r");
    char line[100];
    int i = 0;
    while (fgets(line, 100, file)) {
        length[i] = strlen(line)-1;
        i++;
    }
    fclose(file);
    // write the length[] to a file
    FILE *file2 = fopen("length.txt", "w");
    for (int i = 0; i < 162; i++) {
        fprintf(file2, "%d,\n", length[i]);
    }
    fclose(file2);
    
    return 0;
}