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
void connect_DATA_NODE(DATA_NODE **prev, DATA_NODE **next, DATA_NODE **head) {
    if (*prev) {(*prev) -> next = *next;} // if *prev is not null, then proceed as usual
    else {free(*head); *head = *next;} // if prev is null, then that means next is the head
    *prev = *next;
}


void pack_DATA_NODE(DATA_NODE *node, int8_t zeros, int16_t VAL)
{ 
    uint8_t isNeg = 0, minBits; 
    if (VAL < 0) {VAL = -VAL;  isNeg = 1;}
    minBits = min_bits(VAL);
    node -> rrrrssss = zeros; // number of previous zeros
    node -> rrrrssss <<= 4; // shift to the left 4 bits
    node -> rrrrssss |= minBits; // add on the other side the minimum bits to represent VAL - 1 (even though we know
    // we have eliminated the first bit, so bits = 1 -> 1, bits = 0 -> 0))
    node -> VAL = VAL; // add value on the right of VAL 
    node -> VAL <<= (15 - minBits) + 1; // shift VAL to the left as far as we can, leaving one zero for sign (the +1 is since we dont need leading 1)
    //This means that the maximum value of VAL is +-32767 !! TOTAL_BITSIZE += min_bits_abs(val);
    node -> VAL |= (isNeg << 15); // add the sign
}

