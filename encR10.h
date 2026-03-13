#ifndef __encR10_h__
#define __encR10_h__

#include "globalFile.h"
#include "basisFind.h"
#include "myIO.h"
#include "simulate.h"
#include "decBEC.h"


//global parameters generated and shared by encR10//头文件的作用就是给外部提供接口使用。
extern int K;
extern int S;
extern int H;
extern int n;
extern unsigned int se;
//extern int l;
//extern int L3;
//extern int V;
//extern int a1;
//extern int b1;
//extern int Q;
//extern int A;
//extern int B;
//extern int Y;
//extern int d;
extern VVI inte;//中间符号矩阵
extern VD distri_LT;//cumulative weight distribution for LT symbols
extern VVI R10_LH_adj;
extern VVI R10_LH_val;

//extern double Loss[];
//extern int p;

//number of information bits (excluding seed) contained in an oligo
//equivalently, number of bits in a source symbol
int R10_info_length();

void fun();
void creat_inte(const VVI& info);
VI R10_encode(const VVI &info);
VI R10_xor_VN(unsigned int seed);
#endif
