
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
int max_val(int nBits)
{
    return (1 << nBits) - 1;
}

typedef struct HUFFMAN_NODE
{
    u_int8_t bits;
    HUFFMAN_NODE *left, *right, *parent;
    int weight;

} HUFFMAN_NODE;
/* for AC codes bits = 8 bits, for DC codes bits = 4 bits 
1- assign weights to DATA_NODEs, and store them in HUFFMN_NODEs
2- build the HUFFMAN tree
3- traverse the HUFFMAN and create table of codes
4- encode the DATA_NODEs with the table of codes
*/
typedef struct DATA_NODE 
{
    u_int8_t bits;
    int weight;
} DATA_NODE;

int main () {

    printf("%d\n", max_val(15));
    return 0;
}
 