#ifndef __decBEC_h__
#define __decBEC_h__

#include "globalFile.h"
#include "myIO.h"
#include "triangulation.h"

//global parameters shared by different programs
extern int	m;//#CNs
extern int	n;//#VNs
extern VVI	CN_adj;//neighbors of CNs
extern VVI	VN_adj;//neighbors of VNs
extern VVI	CN_val;//XOR value of CNs

//global parameters generated and shared by decBEC
extern bool use_BP;//true: BP; false: SGE
extern VVI	rec_info;//recovered info: empty for non-recovered

//input: m, n, CN_adj, CN_val (does not change)
//[rec_info] records the recoved info
//return true/false: rec_info has n/<n recovered source symbols
bool decBEC();

#endif
