#define DEBUG true

#include "PDCCH_encoder.h"

uint8_t PDCCH_enc(uint8_t *info, //original data
		  uint32_t info_length,//data length
		  uint32_t encoded_bits,//number of encoded bits, used for rate-matching
		  uint8_t RNTI[16],
		  uint32_t *coded_bits // output of coded bits
		  )
{
    uint8_t  crc_init[25]={1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 1};
    uint8_t P=24; //crc length
    uint32_t N=get_3GPP_N(info_length,encoded_bits,9);

    uint8_t *crc_interleaver_pattern=(uint8_t *)malloc((info_length+P)*sizeof(uint8_t));
    if(get_3GPP_crc_interleaver_pattern(info_length+P,crc_interleaver_pattern)) {
      printf("CRC interleavering error!\n");
      return 1;
    }

    uint32_t *rate_matching_pattern_t=(uint32_t *)malloc(encoded_bits*sizeof(uint32_t));
    mode mode_ibp;
        printf("K=%d,E=%d,N=%d\n",info_length,encoded_bits,N);
    if(get_3GPP_rate_matching_pattern(info_length+P,N,encoded_bits,rate_matching_pattern_t,&mode_ibp)) {
      printf("rate-matching pattern error!\n");
      return 1;
    }
    //    printf("K=%d,E=%d,mode=%d\n",K,E,mode_ibp);

    uint32_t *Q_N=(uint32_t*)malloc(N*sizeof(uint32_t));
    if(get_3GPP_sequence_pattern(Q_N,N)) {
      printf("get sequence pattern error!\n");
      return 1;
    }

    uint32_t *info_bit_pattern_t=(uint32_t*)malloc(N*sizeof(uint32_t));
    if (get_3GPP_info_bit_pattern(info_length+P,Q_N,N,rate_matching_pattern_t,encoded_bits,mode_ibp,info_bit_pattern_t)) {
      printf("get info bit pattern error!\n");
      return 1;
    }
    
    //    uint32_t *e=malloc(E*sizeof(uint32_t));
    if(DS1CA_polar_encoder(info,info_length,crc_init,P,RNTI,crc_interleaver_pattern,info_length,info_bit_pattern_t,N,rate_matching_pattern_t,encoded_bits,coded_bits)) {
      printf("kernal encoder error!\n");
      return 1;
    }
    return 0;
}
