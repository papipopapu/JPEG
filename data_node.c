#include "image_compression.h"



DATA_NODE *new_DATA_NODE()
{
    DATA_NODE* fetus = NULL;
    fetus = (DATA_NODE*)malloc(sizeof(DATA_NODE));
    fetus -> next = NULL;
    return fetus;
}

void free_DATA_NODE_list(DATA_NODE* head)
{
    DATA_NODE* temp;
    while (head != NULL)
    {
        temp = head;
        head = head -> next;
        free(temp);
    }
}
void pack_DATA_NODE(DATA_NODE *node, int8_t zeros, int16_t VAL, size_t *TOTAL_BITSIZE)
{ 
    u_int8_t isNeg = 0, minBits = min_bits_abs(VAL);
    if (VAL < 0) {VAL = -VAL;  isNeg = 1;}
    node -> zeros_bitsVAL = zeros; // number of previous zeros
    node -> zeros_bitsVAL <<= 4; // shift to the left 4 bits
    node -> zeros_bitsVAL |= minBits; // add on the other side the minimum bits to represent VAL - 1 (even though we know
    // we have eliminated the first bit, so bits = 1 -> 1, bits = 0 -> 0))
    node -> VAL = VAL; // add value on the right of VAL 
    node -> VAL <<= (15 - minBits) + 1; // shift VAL to the left as far as we can, leaving one zero for sign (the +1 is since we dont need leading 1)
    //This means that the maximum value of VAL is +-32767 !! TOTAL_BITSIZE += min_bits_abs(val);
    node -> VAL |= (isNeg << 15); // add the sign
    *TOTAL_BITSIZE += minBits + 8;
    // if val = 0, the bits VAL doesnt matter bcs its determined by size = 0
    //so we can just remove the most significant bit, since its always going to be 1

    // possibilities

    /* bit_size = 0 */
    // (0000 0000) (.) EOB, or  -> unlucky, 0 in the last position
    // (1111 0000) (.) 16 chain of zeros
    // (0010 0000) (.) zeros preceded by 3 (n) zeros -> unlucky, 0 in the last position
    /* bit_size != 0 *-> we now take {1 sign bit + bit_size - 1 removed leading bit} = {bit_size} */
    // (0010 0001) (1(1)) -1, (removed leading bit in parenthesis)
    // (0010 0001) (0(1)) 1
    // (0010 0010) (1(1)0) -2
    // (0010 0010) (0(1)1) 3
    // ... etc

    // now we would huffman encode RRRRSSSS


}

