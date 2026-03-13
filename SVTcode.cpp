#include "SVTcode.h"

bool SVT_valid_oligo(const VB& oligo) {
	if (oligo.size() != 2 * N) {
		return	false;
	}
	
	VVB tmp = de_interleave(oligo);
	return	syndrome(tmp[0], 2 * N) == 0 && max_run_binary(tmp[0]) <= P && \
	syndrome(tmp[1], P) == 0 && Hamming_weight(tmp[1]) % 2 == 0;
}

int SVT_info_length() {
	return	2 * N - (int)(log2(N) - eps + 2) - (int)(log2(P) - eps + 2);
}

VB SVT_encode_binary(const VB& info) {
	assert(P > 0 && P <= N);
	assert(N == info.size() + (int)(log2(P) - eps + 2));
	
	int idx = -1;
	int syn = 0;
	VB res(N, false);
	for (int i = 1; i <= N; ++i) {
		if (i > P || (i < P && (i & -i) != i)) {//not a power of 2
			res[i - 1] = info[++idx];
			syn += i * res[i - 1];
		}
	}
	
	syn = (P - syn % P) % P;
	for (int i = 0; (syn >> i) > 0; ++i) {//make syn 0
		if ((syn >> i) & 1) {
			res[(1 << i) - 1] = true;
			syn -= (1 << i);
		}
	}
	
	if (Hamming_weight(res) % 2 != 0) {//make weight even
		res[P - 1] = true;
	}
	
	//	output(info, "SVT_info");
	//	output(res, "SVT_code");
	return	res;
}

VB SVT_decode_binary(const VB& noiseCode, int left, int right) {
	VB res;
	bool b = Hamming_weight(noiseCode) % 2;
	int syn = syndrome(noiseCode, P);
	
	//	cerr << "left right  " << left << " " << right << endl;
	
	//correct single insertion that may happen at the i-th bit of noiseCode, i \in [left, right]
	if (noiseCode.size() == N + 1) {
		if (right - left > P) {//maximum run of original codeword is at most P
			return	res;
		}
		
		int weight = 0;
		for (int i = noiseCode.size() - 1; i >= left; --i) {
			if (i <= right && noiseCode[i] == b && (b * (i + 1) + weight - syn) % P == 0) {
				res = noiseCode;
				res.erase(res.begin() + i);
				return	res;
			}
			weight += noiseCode[i];
		}
	}
	
	//correct single deletion that may happen just before the i-th bit of noiseCode, i \in [left, right]
	if (noiseCode.size() == N - 1) {
		if (right - left + 1 > P) {
			return	res;
		}
		
		int weight = 0;
		for (int i = N - 1; i >= left; --i) {
			if ((syn + weight + b * (i + 1)) % P == 0) {
				res = noiseCode;
				res.insert(res.begin() + i, b);
				return	res;
			}
			if (i > 0) {
				weight += noiseCode[i - 1];
			}
		}
	}
	
	//correct single substitution that may happen at the i-th bit of noiseCode, i \in [left, right]
	if (noiseCode.size() == N) {
		if (syn == 0 && b == false) {
			return	noiseCode;
		}
		
		if (left == right && b == true && (syn + (1 - 2 * noiseCode[left]) * (left + 1)) % P == 0) {
			res = noiseCode;
			res[left] = !res[left];
			return	res;
		}
	}
	
	//failed at this point
	return	VB(0);
}

VB SVT_extract_info_binary(const VB& code) {
	assert(N == code.size());
	
	VB res;
	for (int i = 1; i <= N; ++i) {
		if (i > P || (i < P && (i & -i) != i)) {//not a power of 2
			res.push_back(code[i - 1]);
		}
	}
	return	res;
}

VB SVT_encode(const VB& info) {
	assert(info.size() == SVT_info_length());
	int k = N - (int)(log2(N) - eps + 2);
	VB codeA = L_encode_binary(VB(info.begin(), info.begin() + k));//L encode first part
	if (max_run_binary(codeA) > P) {//the run-length is invalid for SVT code
		return	VB(0);
	}
	VB codeB = SVT_encode_binary(VB(info.begin() + k, info.end()));//SVT encode second part
	return	interleave(codeA, codeB);
}

VB SVT_decode(const VB& noiseOligo) {
	VVB tmp = de_interleave(noiseOligo);
	
	VB codeA = L_decode_binary(tmp[0]);
	if (codeA.empty()) {
		return	VB(0);
	}
	
	PII pos = error_region(codeA, tmp[0]);
	VB codeB = SVT_decode_binary(tmp[1], pos.first, pos.second);
	if (codeB.empty()) {
		return	VB(0);
	}
	
	return	interleave(codeA, codeB);
}

VB SVT_extract_info(const VB& oligo) {
	assert(2 * N == oligo.size());
	
	VVB tmp = de_interleave(oligo);
	VB res = L_extract_info_binary(tmp[0]);
	tmp[1] = SVT_extract_info_binary(tmp[1]);
	res.insert(res.end(), tmp[1].begin(), tmp[1].end());
	return	res;
}

VVI SVT_read(const VVI& write_oligo) {
	assert(write_size == write_oligo.size());
	VI cnt(write_size, 0);
	for (int i = 0; i < read_size; ++i) {
		++cnt[randi(write_size)];//generating sequencing number for each write_oligo
	}
	
	VVI res;
	VVI err;
	for (int i = 0; i < write_size; ++i) {
		VB oligo = decompress(write_oligo[i], 2 * N);
		VI info = compress_ES(SVT_extract_info(oligo));
		
		for (int j = cnt[i]; j > 0; --j) {
			--cnt[i];
			
			VB noiseOligo = add_noise(oligo);
			if (detect_only) {
				if (!SVT_valid_oligo(noiseOligo)) {
					continue;
				}
			}
			else {
				noiseOligo = SVT_decode(noiseOligo);
				if (noiseOligo.empty()) {
					continue;
				}
			}
			
			//noiseOligo is a valid oligo at this point
			VI outInfo = compress_ES(SVT_extract_info(noiseOligo));
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


