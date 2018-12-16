#include "utils.h"

uint8_t PDCCH_enc(uint8_t *info, //original data
		  uint32_t info_length,//data length
		  uint32_t encoded_bits,//number of encoded bits, used for rate-matching
		  uint8_t RNTI[16],
		  uint32_t *coded_bits // output of coded bits
		  );
