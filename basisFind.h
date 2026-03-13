#ifndef BASISFIND_H
#define BASISFIND_H
#include "globalFile.h"
#include "triangulation.h"
#include"encR10.h"

//global parameters shared by different programs
extern int	m;//#CNs
extern int	n;//#VNs
extern int	l;//int-width of CN value
extern VVI	CN_adj;//neighbors of CNs
extern VVI	CN_val;//XOR value of CNs
extern string dec_type;//BP, tri, str, sort <-> BP, triangulation-based BFA, straightforward BFA, sorted-weight BFA

extern  VI  CN_con;//如果是约束符号则为1，否则为0；
//extern  double CON_basis_P;//约束行在基内的平均概率；
//extern  double CON_basis_AP;//所有约束行都在基向量里的概率；
//extern  double B_weight;
extern bool cmp1(int CN_idx1, int CN_idx2);
extern bool cmp(int CN_idx1, int CN_idx2);
extern VVI read_EES;
extern VVI copy_EES;

//global parameters generated and shared by basisFind
extern VB	CN_err;//1: CN is incorrect; 0: CN is correct
extern int	times_n_cor_basis;//#times for having n correct basis elements
//Following parameters for the cases of having n correct basis elements
extern LL	num_inc_basis;//average number of incorrect basis elements 
extern double	basis_wt;//mean of basis elements' weights
extern VD	p_att_mean;//0/1 mean of probability that a correct/incorrect basis element attends LR
extern VD	p_att_var;//0/1 variance of probability that a correct/incorrect basis element attends LR

//input: m, n, l, CN_adj, CN_val
//return 0: n correct basis elements and they uniquely attend top n linear representations; 
//1: error: less than n correct basis elements;
//2: error: n correct basis elements but they do not uniquely attend top n linear representations
int	BFA();

#endif
