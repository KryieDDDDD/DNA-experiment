#include "RScode.h"

int	q = -1;//= (1 << deg) - 1; this variable is used for judging initialization
int	num_rb;//(dmin - 1) * deg: the number of redundant bits
VI	vecExp;//vector: [0...01, 1...11] -> exponent of primitive root: [0, q-1]
VI	expVec;//exponent of primitive root: [0, 2q-1] -> vector: [0...01, 1...11]

//parameters for RS code
VI	gx;//generator polynomial
VI	remErr;//remainder to error: (position << deg) + value

int GF_multiply(int a, int b) {
	assert(a >= 0 && a <= q);
	assert(b >= 0 && b <= q);
	if (a == 0 || b == 0) {
		return	0;
	}
	else {
		return	expVec[vecExp[a] + vecExp[b]];
	}
}

int GF_inverse(int b) {
	assert(b > 0 && b <= q);
	return	expVec[q - vecExp[b]];
}

int GF_divide(int a, int b) {
	return	GF_multiply(a, GF_inverse(b));
}

VI poly_multiply(const VI& A, const VI& B) {
	VI res(A.size() + B.size() - 1, 0);
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < B.size(); ++j) {
			res[i + j] ^= GF_multiply(A[i], B[j]);
		}
	}
	return	res;
}

VI poly_modulo(VI A, const VI& mod) {
	assert(A.size() >= mod.size() && mod.size() > 1 && mod[0] > 0);
	for (int i = 0; i + mod.size() <= A.size(); ++i) {
		if (A[i] != 0) {
			int coe = GF_divide(A[i], mod[0]);
			for (int j = 1; j < mod.size(); ++j) {
				A[i + j] ^= GF_multiply(mod[j], coe);
			}
		}
	}
	return	VI(A.end() - mod.size() + 1, A.end());
}

void RS_init() {
	assert(deg > 0 && deg <= GF_MAX_DEG);
	q = (1 << deg) - 1;//q > 0, initialized
	num_rb = (dmin - 1) * deg;
	assert(dmin > 0 && 2 * N > num_rb && 2 * N <= q * deg);//RS code matches oligo
	
	vecExp.resize(q + 1);
	expVec.resize(2 * q);
	vecExp[0] = -1;//vec: 0 -> exp: -1 (undefined)
	int vec = 1;//primitive_root^0 <-> 00...01
	for (int i = 0; i < expVec.size(); ++i) {
		if (i < q) {
			vecExp[vec] = i;
		}
		expVec[i] = vec;
		vec <<= 1;
		if (vec >> deg) {
			vec ^= primPoly[deg];
		}
	}
	
	gx.resize(1, 1);
	for (int i = 1; i < dmin; ++i) {
		VI factor(2, 1);
		factor[1] = expVec[i];//x - primitive_root^i
		gx = poly_multiply(gx, factor);
	}
	
	if (dmin != 3 || num_rb > 20) {//cannot generate remErr to correct errors
		return;
	}
	remErr.resize(1 << num_rb, -1);
	remErr[0] = 0;
	VI code((2 * N + deg - 1) / deg, 0);
	for (int i = 0; i < code.size(); ++i) {
		for (int v = 1; v <= q; ++v) {
			code[i] = v;
			VI rem = poly_modulo(code, gx);
			int idx = 0;
			for (int j = 0; j < rem.size(); ++j) {
				idx = (idx << deg) + rem[j];
			}
			remErr[idx] = (i << deg) + v;
		}
		code[i] = 0;
	}
}

bool RS_valid(const VB& oligo) {
	if (q < 0) {
		RS_init();
	}
	
	if (oligo.size() != 2 * N) {
		return	false;
	}
	
	VB info = VB(oligo.begin(), oligo.end() - num_rb);
	return	RS_encode(info) == oligo;
}

int RS_info_length() {
	if (q < 0) {
		RS_init();
	}
	
	return	2 * N - num_rb;
}

VB RS_encode(const VB& info) {
	if (q < 0) {
		RS_init();
	}
	
	assert(info.size() == RS_info_length());
	
	VI code = compress(info, deg);//(info.size() % deg) 0s are appended, but will not be returned
	code.resize(code.size() + dmin - 1, 0);//append (dmin - 1) 0s
	VI rem = poly_modulo(code, gx);
	VB rem_bit = decompress(rem, rem.size() * deg, deg);
	VB res = info;
	res.insert(res.end(), rem_bit.begin(), rem_bit.end());
	return	res;
}

VB RS_decode(const VB& noiseCode) {
	if (q < 0) {
		RS_init();
	}
	
	assert(!remErr.empty());
	
	if (noiseCode.size() != 2 * N) {
		return	VB(0);
	}
	
	VI code = compress(VB(noiseCode.begin(), noiseCode.end() - num_rb), deg);
	VI rem = compress(VB(noiseCode.end() - num_rb, noiseCode.end()), deg);
	code.insert(code.end(), rem.begin(), rem.end());
	
	rem = poly_modulo(code, gx);
	int idx = 0;
	for (int j = 0; j < rem.size(); ++j) {
		idx = (idx << deg) + rem[j];
	}
	if (remErr[idx] < 0) {
		return	VB(0);
	}
	else {
		code[remErr[idx] >> deg] ^= (remErr[idx] & q);
		VB res = decompress(code, code.size() * deg, deg);
		res.erase(res.begin() + 2 * N - num_rb, res.end() - num_rb);
		return	res;
	}
}

VB RS_extract_info(const VB& code) {
	return	VB(code.begin(), code.begin() + RS_info_length());
}

VVI RS_read(const VVI& write_oligo) {
	if (q < 0) {
		RS_init();
	}
	
	assert(write_size == write_oligo.size());
	VI cnt(write_size, 0);
	for (int i = 0; i < read_size; ++i) {
		++cnt[randi(write_size)];//generating sequencing number for each write_oligo
	}
	
	VVI res;
	VVI err;
	for (int i = 0; i < write_size; ++i) {
		VB oligo = decompress(write_oligo[i], 2 * N);
		VI info = compress_ES(RS_extract_info(oligo));
		
		for (int j = cnt[i]; j > 0; --j) {//如果这里cnt[i]=0;就不会进入循环，于是oligo被擦除；
			--cnt[i];
			
			VB noiseOligo = add_noise(oligo);
			if (detect_only) {
				if (!RS_valid(noiseOligo)) {
					continue;
				}
			}
			else {
				noiseOligo = RS_decode(noiseOligo);
				if (noiseOligo.empty()) {
					continue;
				}
			}
			
			//noiseOligo is a valid oligo at this point
			VI outInfo = compress_ES(RS_extract_info(noiseOligo));
			if (outInfo == info) {
				++cnt[i];
			}
			else {//ignore the situation of becoming a correct encoded symbol
				err.push_back(outInfo);
			}
		}
		
		if (cnt[i] > 0) {
			res.push_back(info);
			res.back().push_back(cnt[i]);
		}
	}
	
	sort(err.begin(), err.end());
	for (int i = 0; i < err.size();) {
		int j = i;
		while (++j < err.size() && err[j] == err[i]);
		err[i].push_back(i - j);//negative for incorrect encoded symbols
		res.push_back(err[i]);
		i = j;
	}
	
	return	res;
}
