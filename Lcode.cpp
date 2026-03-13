#include "Lcode.h"

VS originalStrand;

int syndrome(const VB& A, int mod) {
	int res = 0;
	for (int i = 0; i < A.size(); ++i) {
		res += A[i] * (i + 1);
	}
	return	res % mod;
}

PII error_region(const VB& code, const VB& noiseCode) {
	int len = min(code.size(), noiseCode.size());
	//leftmost and rightmost possible position noiseCode has an error: 
	//insert [0, N]; delete [0, N-1]; substitution [0, N-1] ; correct <-1, N>
	int left, right;
	for (left = 1; left <= len; ++left) {
		if (code[code.size() - left] != noiseCode[noiseCode.size() - left]) {
			break;
		}
	}
	left = code.size() - left + (noiseCode.size() > code.size());
	for (right = 0; right < len; ++right) {
		if (code[right] != noiseCode[right]) {
			break;
		}
	}
	return	make_pair(left, right);
}

bool L_valid_oligo(const VB& oligo) {
	if (oligo.size() != 2 * N) {
		return	false;
	}
	
	VVB tmp = de_interleave(oligo);
	return	syndrome(tmp[0], 2 * N) == 0 && syndrome(tmp[1], 2 * N) == 0;
}

int L_info_length() {
	return	2 * (N - (int)(log2(N) - eps + 2));
}

VB L_encode_binary(const VB& info) {
	assert(N == info.size() + (int)(log2(N) - eps + 2));
	
	int idx = -1;
	int syn = 0;
	int mod = 2 * N;
	VB res(N, false);
	for (int i = 1; i < N; ++i) {
		if ((i & -i) != i) {//not a power of 2
			res[i - 1] = info[++idx];
			syn += i * res[i - 1];
		}
	}
	
	syn = (mod - syn % mod) % mod;
	if (syn >= N) {
		res[N - 1] = true;
		syn -= N;
	}
	for (int i = 0; (syn >> i) > 0; ++i) {
		if ((syn >> i) & 1) {
			res[(1 << i) - 1] = true;
			syn -= (1 << i);
		}
	}
	
	//	output(info, "L_info");
	//	output(res, "L_code");
	return	res;
}

VB L_decode_binary(const VB& noiseCode) {
	VB res;
	int mod = 2 * N;
	int syn = syndrome(noiseCode, mod);
	
	//correct single insertion
	if (noiseCode.size() == N + 1) {
		int weight = 0;
		for (int i = noiseCode.size(); i > 0; --i) {
			if (i * noiseCode[i - 1] + weight == syn) {
				res = noiseCode;
				res.erase(res.begin() + i - 1);
				return	res;
			}
			weight += noiseCode[i - 1];
		}
	}
	
	//correct single deletion
	if (noiseCode.size() == N - 1) {
		int weight = 0;
		for (int i = N; i > 0; --i) {
			if ((syn + weight) % mod == 0 || (syn + weight + i) % mod == 0) {
				res = noiseCode;
				res.insert(res.begin() + i - 1, (syn + weight + i) % mod == 0);
				return	res;
			}
			if (i > 1) {
				weight += noiseCode[i - 2];
			}
		}
	}
	
	//correct single substitution
	if (noiseCode.size() == N) {
		if (syn == 0) {
			return	noiseCode;
		}
		
		if (syn <= N && noiseCode[syn - 1] == true) {//flip the (syn-1)-th bit from 1 to 0
			res = noiseCode;
			res[syn - 1] = false;
			return	res;
		}
		
		if (syn >= N && noiseCode[mod - syn - 1] == false) {//flip the (mod-syn-1)-th bit from 0 to 1
			res = noiseCode;
			res[mod - syn - 1] = true;
			return	res;
		}
	}
	
	//failed at this point
	return	VB(0);
}

VB L_extract_info_binary(const VB& code) {
	assert(N == code.size());
	
	VB res;
	for (int i = 1; i < N; ++i) {
		if ((i & -i) != i) {//not a power of 2
			res.push_back(code[i - 1]);
		}
	}
	return	res;
}

VB L_encode(const VB& info) {
	int k = L_info_length();
	assert(info.size() == k);
	VVB info2 = de_interleave(info);
	VB codeA = L_encode_binary(info2[0]);//encode 0, 2, 4th ... bits
	VB codeB = L_encode_binary(info2[1]);//encode 1, 3, 5th ... bits
	return	interleave(codeA, codeB);
}

VB L_decode(const VB& noiseOligo) {
	VVB tmp = de_interleave(noiseOligo);
	VVB code(tmp.size());
	int left = -1, right = N;
	for (int i = 0; i < tmp.size(); ++i) {
		code[i] = L_decode_binary(tmp[i]);
		if (code[i].empty()) {
			return	code[i];
		}
		
		PII pr = error_region(code[i], tmp[i]);
		left = max(left, pr.first);
		right = min(right, pr.second);
	}
	
	if (left <= right) {
		return	interleave(code[0], code[1]);
	}
	else {
		return	VB(0);
	}
}

VB L_extract_info(const VB& oligo) {
	assert(2 * N == oligo.size());
	
	VVB tmp = de_interleave(oligo);
	VB res = L_extract_info_binary(tmp[0]);
	tmp[1] = L_extract_info_binary(tmp[1]);
	//res.insert(res.end(), tmp[1].begin(), tmp[1].end());
	res = interleave(res, tmp[1]);
	return	res;
}

VVI L_read(const VVI& write_oligo) {
	assert(write_size == write_oligo.size());
	VI cnt(write_size, 0);
	for (int i = 0; i < read_size; ++i) {
		++cnt[randi(write_size)];//generating sequencing number for each write_oligo
	}
	
	VVI res;
	VVI err;
	for (int i = 0; i < write_size; ++i) {
		VB oligo = decompress(write_oligo[i], 2 * N);
		VI info = compress_ES(L_extract_info(oligo));
		
		for (int j = cnt[i]; j > 0; --j) {
			--cnt[i];
			
			VB noiseOligo = add_noise(oligo);
			if (detect_only) {
				if (!L_valid_oligo(noiseOligo)) {
					continue;
				}
			}
			else {
				noiseOligo = L_decode(noiseOligo);
				if (noiseOligo.empty()) {
					continue;
				}
			}
			
			//noiseOligo is a valid oligo at this point
			VI outInfo = compress_ES(L_extract_info(noiseOligo));
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

int seed_cmp(string s1, string s2) {//return 1 if s1 > s2 (A < T < C < G); 0 if =, -1 if <
	//assert(s1.size() == 21 && s2.size() == 21);
	int i = 0;
	for (i = 0; i < s1.size(); ++i)
		if (s1[i] != s2[i])
			break;
	if (i == s1.size())
		return 0;

	if (s1[i] == 'G')
		return 1;
	else if (s1[i] == 'C') {
		if (s2[i] == 'G')
			return -1;
		else return 1;
	}
	else if (s1[i] == 'T') {
		if (s2[i] == 'A')
			return 1;
		else
			return -1;
	}
	else
		return -1;

}

int bin_search(string s1) {
	int left = 0;
	int right = originalStrand.size() - 1;
	string seed1 = string(s1.begin(), s1.begin() + 1) + string(s1.begin() + 3, s1.begin() + 4) + string(s1.begin() + 5, s1.begin() + 8) + string(s1.begin() + 9, s1.begin() + 13);
	while (left <= right) {
		int mid = (left + right) / 2;
		//string seed2 = string(seedInfo[mid].first.begin(), seedInfo[mid].first.begin() + 21);
		string stemp = originalStrand[mid];
		string seed2 = string(stemp.begin(), stemp.begin() + 1) + string(stemp.begin() + 3, stemp.begin() + 4) + string(stemp.begin() + 5, stemp.begin() + 8) + string(stemp.begin() + 9, stemp.begin() + 13);
		
		int itemp = seed_cmp(seed2, seed1);
		if (itemp == 0) {
			if (s1 == originalStrand[mid])
				return mid;
			else
				return -1;
		}
		else if (itemp > 0)
			right = mid - 1;
		else
			left = mid + 1;
	}
	return -1;
}



void L_decode_from_txt()
{
	ifstream ifs1, ifs2;
	ofstream ofs1, ofs2;
	//ifs1.open("C:\\Users\\xinmatrix\\Desktop\\Data\\1.txt", ios::in);
	//ifs2.open("C:\\Users\\xinmatrix\\Desktop\\Data\\C1.txt", ios::in);
	ifs1.open("C:\\Users\\xinmatrix\\Desktop\\Data\\1.txt", ios::in);
	ifs2.open("C:\\Users\\xinmatrix\\Desktop\\Data\\C1.txt", ios::in);
	ofs1.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\1_detect_only.txt", ios::out);
	ofs2.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\1_detect_info.txt", ios::out);

	//if (!ifs1.is_open() || !ifs2.is_open() || !ofs1.is_open() || !ofs2.is_open())
		//assert(0);
	int t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;//t1 - no error; t2 - dec suc; t3 - undetec error; t4 - no dec res; t5 - < 249 || > 251
	int a1 = 0, a2 = 0;
	long long dis = 0;
	string stemp;
	char ctemp;
	int maxLine = 12173158;//0 - 11702066; 1 - 12173158
	time_t startTime = clock();
	while (getline(ifs2, stemp))
		originalStrand.push_back(stemp);
	int lineNum = 0, decNum = 0;
	while (getline(ifs1, stemp)) {
		++lineNum;
		if (!(lineNum % 100000)) {
			cout << 1.0 * lineNum / maxLine << ", " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
			
		}
		//if (stemp.size() < N + 1 || stemp.size() > N + 3) {
		if (stemp.size() != N + 2) {/////////////////////////////////////////////////for detect only
			++t5;
			continue;
		}

		++decNum;
		ctemp = stemp[0];
		//cout << stemp;
		stemp = string(stemp.begin() + 1, stemp.end() - 1);
		VB stempBin = string2bin(stemp);
		VB decBin = L_decode(stempBin);
		bool flag_same = (decBin == stempBin) ? 1 : 0;
		if (!flag_same) {////////////////////////////////////////////////////////////for detect only 
			++a1;
			continue;
			
		}
		if (decBin.empty()) {
			++t4;
			continue;
		}
		
		stemp = bin2string(decBin);
		stemp = ctemp + stemp + ctemp;
		//cout << stemp;
		//string seed = string(stemp.begin(), stemp.begin() + 1) + string(stemp.begin() + 3, stemp.begin() + 4) + string(stemp.begin() + 5, stemp.begin() + 8) + string(stemp.begin() + 9, stemp.begin() + 13);
		int itemp = bin_search(stemp);
		if (itemp < 0) {
			++t3;
			//cout << "t3++, " << lineNum << endl;
		}
		else {
			if (flag_same)
				++t1;
			else
				++t2;
		}
		decBin = L_extract_info(decBin);
		VB seed = VB(decBin.begin(), decBin.begin() + bw_seed);
		reverse(seed.begin(), seed.end());
		decBin = VB(decBin.begin() + bw_seed, decBin.end());
		decBin = SRT_decode(decBin);
		if (decBin.empty()) {
			++t6;
			--t3;
			//cout << "t6++, " << lineNum << endl;
			continue;
		}

		decBin.insert(decBin.begin(), seed.begin(), seed.end());
		VI dec = compress_ES(decBin);
		++a2;
		ofs1 << (itemp < 0) << " ";
		for (int i = 0; i < dec.size(); ++i)
			ofs1 << dec[i] << " ";
		ofs1 << endl;
		
	}
	ofs2 << "process strands : " << lineNum << endl;
	ofs2 << "L distance > 1 : " << t5 << endl;
	ofs2 << "received with no error : " << t1 << endl;
	ofs2 << "SRT dec error : " << t6 << endl;
	ofs2 << "have errors & dec successfully : " << t2 << endl;
	ofs2 << "have errors & undet error " << t3 << endl;
	ofs2 << "have errors & no dec result : " << t4 << endl;
	ofs2 << "a1 = " << a1 << endl;
	ofs2 << "a2 = " << a2 << endl;
}


