//Single Parity Check Code

#ifndef __SPCC_h__
#define __SPCC_h__

#include "globalFile.h"
#include "myIO.h"
#include "simulate.h"

VB SPCC_encode(VB const& info);

int SPCC_info_length();

bool valid_SPCC_cw(VI& cw);

VVB SPCC_dec_4A_seq(VB cw);

VB SPCC_dec_4A_cluster(const VVB& cluster, int& freq);

VB SPCC_extract_info(VB& cw);

#endif