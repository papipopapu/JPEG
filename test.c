
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
int max_val(int nBits)
{
    return (1 << nBits) - 1;
}
int main () {

printf("%d\n", max_val(15));
    return 0;
}
 