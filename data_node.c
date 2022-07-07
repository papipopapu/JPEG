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
void pack_DATA_NODE(DATA_NODE *node, int8_t zeros, int16_t VAL)
{ 
    u_int8_t isNeg = 0, minBits = min_bits_abs(VAL);
    if (VAL < 0) {VAL = -VAL;  isNeg = 1;}
    node -> zeros_bitsVAL = zeros; // number of previous zeros
    node -> zeros_bitsVAL <<= 4; // shift to the left 4 bits
    node -> zeros_bitsVAL |= minBits; // add on the other side the minimum bits to represent VAL - 1 (a
    //there always is a sign value in the first bit of VAL
    node -> VAL = VAL; // add value on the right of VAL 
    node -> VAL <<= (15 - minBits); // shift VAL to the left as far as we can, leaving one zero for sign 
    //This means that the maximum value of VAL is +-32767 !! 
    node -> VAL |= (isNeg << 15); // add the sign

}