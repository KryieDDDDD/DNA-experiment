#include "encR10.h"

int K;
unsigned int se = 0;
int L3;
int S;
int H;
int H1;
VVI inte;
VD distri_LT;
VVI R10_LH_adj;
VVI R10_LH_val;

//定义组合函数
int Z(int a, int b) {
	double sum = 0;
	for (int i = b + 1; i <= a; ++i) {
		sum += log10(i);
	}
	for (int i = 1; i <= a - b; ++i) {
		sum -= log10(i);
	}
	double result = pow(10, sum);
	return int(result);
}

void fun() {
	int i, j, W;
	//求W
	for (i = 1;; i++) {
		if (i * (i - 1) >= 2 * K) {
			W = i;
			break;
		}
		else {
			continue;
		}
	}
	//	int W = i;
	//求S

	for (i = ceil(0.01 * K) + W;; i++)
	{
		for (j = 2; j < i; j++) {
			if (i % j == 0) {
				break;
			}
		}
		if (i <= j) {
			S = i;
			break;
		}
	}
	//	for(i = ceil(0.01 * K) + W;; i ++){
	//		for(j = 2; j <= (i/j); j++){//一个数如果不能被2到它的平方根的所有数整除，它就是质数
	//			if(!(i%j)){
	//				break;//如果找到，就不是质数
	//			}
	//		}
	//		if(j > (i/j)){
	//			S = i;
	//			break;
	//		}
	//	}
	//求H、n	
	for (i = 1;; i++) {
		if (Z(i, int(ceil(double(i) / 2))) >= K + S) {//注意两整数相除得到的结果是小数的整数部分，故需要在前面加double;
			H = i;
			break;
		}
		else {
			continue;
		}
	}
	//	H = i;
	H1 = ceil(double(i) / 2);
	n = K + S + H;
	//求L'

	for (i = n;; i++) {
		for (j = 2; j < i; j++) {
			if (i % j == 0) {
				break;//如果找到，就不是质数
			}
		}
		if (i <= j) {
			L3 = i;
			break;
		}
	}
}

//由中间符号生成的编码符号的信息位长度
int R10_info_length() {
	if (inner_code == 'R') {
		return	RS_info_length() - bw_seed;
	}
	else if (inner_code == 'L') {
		return	L_info_length() - bw_seed;
	}
	else if (inner_code == 'S') {
		return	SVT_info_length() - bw_seed;
	}
	else if (inner_code == 'P') {
		return SPCC_info_length() - bw_seed;
	}
	else {
		cerr << "inner code is incorrect." << endl;
		assert(false);
	}
}




//求G_LDPC
VVI creat_G_L() {
	VVI G_LDPC(S, VI(K, 0));
	int j = 0;
	for (int i = 0; i < (K + S - 1) / S; i++) {//或者i < ceil(double(K) / S),计算机先计算K/S（整数相除取下整）,再ceil，若不加double,计算值会错误;
		G_LDPC[0][j] = G_LDPC[(i + 1) % S][j] = G_LDPC[2 * (i + 1) % S][j] = 1;
		for (int c = j + 1; c < j + S; c++) {
			if (c < K) {
				for (int r = 0; r < S; r++) {
					G_LDPC[r][c] = G_LDPC[(r + S - 1) % S][c - 1];
				}
			}
			else {
				break;
			}

		}
		j += S;

	}
	return G_LDPC;
}
//求G_Half
VVB creat_G_H() {
	//	fun();
	VVB G_Half(H, VB(K + S, 0));
	int g1_size = 0, m1_size = 0;
	VI m1, g1;
	for (int i = 0;; i++) {
		g1.push_back(i ^ (i / 2));//取整函数返回的是一个double类型的函数，故要转为int型，但是向下取整也可以直接M/N，向上取整可直接用公式（M + N - 1）/N。
		for (int e = 0; e < H && g1[i] != 0; e++) {
			if (g1[i] & (1 << e)) {
				++g1_size;
			}
			else {
				continue;
			}
		}
		if (g1_size == H1) {
			m1.push_back(g1[i]);
			++m1_size;
			g1_size = 0;
		}
		else {
			g1_size = 0;
			continue;
		}
		if (m1_size == K + S) {
			break;
		}
		else {
			continue;
		}

	}
	for (int j = 0; j < K + S; j++) {
		for (int i = 0; i < H; i++) {
			int result;
			result = (m1[j] & (1 << i)) >> i;
			G_Half[i][j] = result;
		}
	}
	return G_Half;
}

//求前S+H行CN_adj
VVI creat_adj() {
	VVI LH_adj;
	LH_adj.resize(S + H);//这里要先规定size，不然无法对LH_adj[i]进行push_back操作，下标越界。
	VVI I_S(S, VI(S, 0));
	VVI I_H(H, VI(H, 0));
	VVB G_Half = creat_G_H();
	VVI G_LDPC = creat_G_L();
	for (int i = 0; i < S; i++) {
		I_S[i][i] = 1;
	}
	for (int i = 0; i < H; i++) {
		I_H[i][i] = 1;
	}
	for (int i = 0; i < S; ++i) {
		for (int j = 0; j < K; ++j) {
			if (G_LDPC[i][j]) {
				LH_adj[i].push_back(j);
			}
		}
		for (int j = K; j < K + S; ++j) {
			if (I_S[i][j - K]) {
				LH_adj[i].push_back(j);
			}
		}
		//		for(int j = K + S; j < K + S + H; ++ j){
		//			if(Z_SH[i][j - K -S]){
		//				CN_adj[i].push_back(j);
		//			}
		//		}
		//		CN_val[i] = inte[CN_adj[i][0]];
		//		for(int j = 1; j < CN_adj[i].size(); ++ j){
		//			//			assert(CN_val[i].size() >= C[CN_adj[i][j]].size());
		//			xorEq(CN_val[i], inte[CN_adj[i][j]]);
		//		}
	}
	for (int i = S; i < S + H; ++i) {
		for (int j = 0; j < K + S; ++j) {
			if (G_Half[i - S][j]) {
				LH_adj[i].push_back(j);
			}
		}
		for (int j = K + S; j < K + S + H; ++j) {
			if (I_H[i - S][j - K - S]) {
				LH_adj[i].push_back(j);
			}
		}
		//		CN_val[i] = inte[CN_adj[i][0]];
		//		for(int j = 1; j < CN_adj[i].size(); ++ j){
		//			//			assert(CN_val[i].size() >= C[CN_adj[i][j]].size());
		//			xorEq(CN_val[i], inte[CN_adj[i][j]]);
		//		}
	}
	return LH_adj;
}


//求前S+H行CN_val
VVI creat_val(VVI& LH_adj) {
	VVI LH_val;
	LH_val.resize(S + H);
	for (int i = 0; i < S + H; ++i) {
		LH_val[i] = inte[LH_adj[i][0]];
		for (int j = 1; j < LH_adj[i].size(); ++j) {
			xorEq(LH_val[i], inte[LH_adj[i][j]]);
		}
	}
	return LH_val;
}

//生成中间符号；
void creat_inte(const VVI& info) {
	inte.clear();
	inte.resize(n);
	for (int i = 0; i < info.size(); i++) {
		inte[i] = info[i];
	}
	for (int i = K; i < K + S + H; ++i) {
		for (int j = 0; j < inte[0].size(); ++j) {
			inte[i].push_back(0);
		}
	}
	//	生成LDPC符号
	for (int i = 0; i < K; i++) {
		int a2 = 1 + ((i / S) % (S - 1));
		int b2 = i % S;
		xorEq(inte[K + b2], inte[i]);//^表示位异或操作，不能直接用在向量之间的异或，&表示按位与。
		b2 = (b2 + a2) % S;
		xorEq(inte[K + b2], inte[i]);
		b2 = (b2 + a2) % S;
		xorEq(inte[K + b2], inte[i]);
	}

	//生成Half符号；
	int g1_size = 0, m1_size = 0;
	VI m1, g1;
	for (int i = 0;; i++) {
		g1.push_back(i ^ (i / 2));
		for (int e = 0; e < H && g1[i] != 0; e++) {
			if (g1[i] & (1 << e)) {
				++g1_size;
			}
			else {
				continue;
			}
		}
		if (g1_size == H1) {
			m1.push_back(g1[i]);
			++m1_size;
			g1_size = 0;
		}
		else {
			g1_size = 0;
			continue;
		}
		if (m1_size == K + S) {
			break;
		}
		else {
			continue;
		}
	}
	for (int h = 0; h < H; h++) {
		for (int j = 0; j < K + S; j++) {
			if (m1[j] & (1 << h)) {
				xorEq(inte[h + K + S], inte[j]);
			}
			else {
				continue;
			}
		}
	}
}

//生成随机种子；

unsigned int creat_seed() {//1 + x^25 + x^26 + x^30 + x^32
	//static unsigned int s = rand32();
	//if (bw_seed == 32) {
	//	s = ((s ^ (s >> 25) ^ (s >> 26) ^ (s >> 30)) << 31) ^ (s >> 1);
	//}
	//else {//< 32
	//	s = (s + 1) % (1u << bw_seed);
	//}
	unsigned int s = 0;
	if (max_run == 3 && inner_code == 'L') {

	}
	else if (max_run == 3 && inner_code == 'P') {
		s = (1u << (bw_seed - 2)) ^ s;
	}
	else if (max_run == 4 && inner_code == 'L') {
		s = (1u << (bw_seed - 1)) ^ s;
	}
	else if (max_run == 4 && inner_code == 'P') {
		s = (1u << (bw_seed - 1)) ^ (1u << (bw_seed - 2)) ^ s;
	}
	return	s;
}

//定义重量分布
int rand_weight() {
	//	binomial_distribution<int> bd(distri_LT.size() - 1, 0.5);
	//	return	bd(PRNG);

	double	pu = randu();
	for (int i = 1; i < distri_LT.size(); ++i) {//as the sum of weights is small, this is faster than binary search
		if (pu <= distri_LT[i]) {
			return	i;
		}
	}
}


//构建CN_adj
VI R10_xor_VN(unsigned int seed) {
	PRNG.seed(seed);
	return	a_choose_b(n, rand_weight());
}
//R10编码得到校验信息与异或信息；



VI R10_encode(const VVI& info) {
	static int seed = -1;

	if (inte.empty()) {
		cout << "seed = " << seed << endl;
		seed = -1;
		creat_inte(info);
		R10_LH_adj = creat_adj();
		R10_LH_val = creat_val(R10_LH_adj);
	}
	while (!valid_run_int(++seed));
	assert(seed < (1 << bw_seed ));

	if (!(seed % 10000))
		cout << "--------" << seed << endl;
	VI VN = R10_xor_VN(seed);
	VI res = inte[VN[0]];
	for (int j = 1; j < VN.size(); ++j) {
		xorEq(res, inte[VN[j]]);
	}
	res.insert(res.begin(), seed);

	return	res;
}