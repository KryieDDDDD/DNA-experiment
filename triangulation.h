#ifndef __triangulation_h__
#define __triangulation_h__

#include "globalFile.h"
#include "myIO.h"

#include "basisFind.h"

//global parameters shared by different programs
extern int	m;//#CNs
extern int	n;//#VNs
extern VVI	CN_adj;//neighbors of CNs
extern VVI	VN_adj;//neighbors of VNs

//global parameters generated and shared by triangulation
extern	bool use_DSDS;//true: use DSDS (max-component); 
					  //false: minimum-weight row -> maximum accumulate column weight
extern	int	num_rec;//number of recovered VNs (size of L)
extern	int	num_inact;//number of inactive VNs (size of R)
extern	VI	CN_pos_idx;//CN pos (stores) idx
extern	VI	VN_idx_pos;//CN idx (is stored at) pos
extern	clock_t	time_triang;//clock time for triangulation

extern int tri_pq;

extern bool use_pq;
extern VD pq_key;

//input: m, n, CN_adj (does not change)
//get form [L R], where L is a lower triangular matrix, |L| = num_rec, |R| = num_inact
//permutation on rows/columns are recorded by CN_pos_idx and VN_idx_pos
//return true: num_rec + num_inact = n; false: otherwise
bool triangulation();
	
#endif
