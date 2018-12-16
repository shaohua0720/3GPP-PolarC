#define DEBUG true
#include "PDCCH_encoder.h"

int main(void)
{
    uint8_t a[15]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    uint8_t RNTI[16]={1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    uint32_t *rs=(uint32_t *)malloc(6*12*2*sizeof(uint32_t));
    if(PDCCH_enc(a,15,6*12*2,RNTI,rs)) {
      printf("encoder error!\n");
    }

    //show the result.
    uint32_t i=0;
    for(i=0;i<6*12*2;i++) {
      printf("%d ",rs[i]);
    }
    return 0;
}
