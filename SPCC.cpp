#include "SPCC.h"

VB SPCC_encode(VB const& info)
{
	assert(info.size() == 2 * (N - 2));
	VI code = bin2qua(info);
	for (int i = 0; i < 2; ++i) {
		code.push_back(0);
		for (int j = code.size() - 3; j >= 0; j -= 2) {
			code.back() ^= code[j];
		}
	}
	return qua2bin(code);
}

int SPCC_info_length() {
	return 2 * N - 4;
}

VVB SPCC_dec_4A_seq(VB cw)
{
	//assert(codewd.size() == 2 * N || codewd.size() == 2 * (N - 1) || codewd.size() == 2 * (N + 1));
	VI cw_int = bin2qua(cw);
	VVB rst;
	vector<pair<int, int>> pos;
	if (cw_int.size() < (N - 1) || cw_int.size() > (N + 1))
		return VVB{};
	else if (cw_int.size() == N) {
		if (valid_SPCC_cw(cw_int)) {
			//cw.resize(cw.size() - 4);
			return VVB{ cw };
		}
		else {
			int odd = 0, even = 0;
			for (int i = 0; i < cw_int.size(); i += 2)
				even ^= cw_int[i];
			for (int i = 1; i < cw_int.size(); i += 2)
				odd ^= cw_int[i];
			if (even && odd)
				return VVB(0);			
			else if (!odd) {
				//cw_int.resize(cw_int.size() - 2);
				for (int i = 0; i < cw_int.size(); i += 2) {
					cw_int[i] ^= even;
					rst.push_back(qua2bin(cw_int));
					cw_int[i] ^= even;
				}
				return rst;
			}
			else {
				//cw_int.resize(cw_int.size() - 2);
				for (int i = 1; i < cw_int.size(); i += 2) {
					cw_int[i] ^= odd;
					rst.push_back(qua2bin(cw_int));
					cw_int[i] ^= odd;
				}
				return rst;
			}
		}
	}
	else if (cw_int.size() == N - 1) {
		for (int i = 0; i < cw_int.size(); ++i) {
			for (int base = 0; base < 4; ++base) {
				if (i > 0 && base == cw_int[i - 1])
					continue;
				VI cw_tmp = cw_int;
				cw_tmp.emplace(cw_tmp.begin() + i, base);
				if (valid_SPCC_cw(cw_tmp)) {
					//cw_tmp.resize(cw_tmp.size() - 2);
					rst.push_back(qua2bin(cw_tmp));
					break;
				}
			}	
		}
		return rst;
	}
	else {
		bool flag = 0;
		for (int i = 0; i < cw_int.size(); ++i) {
			//cout << i << " ";
			if (flag) {
				while (cw_int[i] == cw_int[i - 1]) {
					if (i + 1 < cw_int.size())
						++i;
					else
						return rst;
				}
			}
			VI cw_tmp = cw_int;
			cw_tmp.erase(cw_tmp.begin() + i);
			flag = 0;
			if (valid_SPCC_cw(cw_tmp)) {
				flag = 1;
				//cw_tmp.resize(cw_tmp.size() - 2);
				rst.push_back(qua2bin(cw_tmp));
			}
		}
		return rst;
	}
	
}

VB SPCC_dec_4A_cluster(const VVB& cluster, int& freq) {
	//VVB candidate;
	freq = 0;
	unordered_map<VB, int> candidate;
	for (int i = 0; i < cluster.size(); ++i) {
		VVB list = SPCC_dec_4A_seq(cluster[i]);
		for (auto& c : list) {
			//cout << bin2string(c) << endl;
			auto iter = candidate.find(c);
			if (iter != candidate.end())
				++iter->second;
			else {
				candidate[c] = 1;
				
			}
		}
	}
	vector<pair<VB, int>> candidate_tmp;
	for (const auto& i : candidate)
		candidate_tmp.push_back(i);
	int highest = 0, sec_highest = 0;
	if (candidate_tmp.empty())
		return VB{};
	if (candidate_tmp.size() == 1) {
		freq = candidate_tmp.back().second;
		return candidate_tmp.back().first;
	}
	//if (candidate.size() == 2) {
	//	if (freq[0] == freq[1])
	//		return VB{};
	//	else
	//		return freq[0] > freq[1] ? candidate[0] : candidate[1];
	//}
	highest = candidate_tmp[0].second > candidate_tmp[1].second ? 0 : 1;
	sec_highest = candidate_tmp[0].second < candidate_tmp[1].second ? 0 : 1;
	for (int i = 2; i < candidate_tmp.size(); ++i) {
		if (candidate_tmp[highest].second < candidate_tmp[i].second) {
			sec_highest = highest;
			highest = i;
		}
		else if (candidate_tmp[sec_highest].second < candidate_tmp[i].second)
			sec_highest = i;
	}
	if (candidate_tmp[highest].second == candidate_tmp[sec_highest].second)
		return VB{};
	freq = candidate_tmp[highest].second;
	return candidate_tmp[highest].first;
}

VB SPCC_extract_info(VB& cw)
{
	cw.resize(cw.size() - 4);
	return cw;
}

bool valid_SPCC_cw(VI& cw)
{
	//assert(cw.size() == N);
	if (cw.size() != N)
		return false;
		//cerr << "SPCC.cpp -> valid_SPCC_cw input length != N" << endl;
	VI cw_tmp = cw;
	//VI cw = bin2qua(codewd);
	for (int i = cw_tmp.size() - 3; i >= 0; i -= 2) {
		cw_tmp.back() ^= cw_tmp[i];
	}
	if (cw_tmp.back() != 0)
		return false;
	cw_tmp.pop_back();
	for (int i = cw_tmp.size() - 3; i >= 0; i -= 2) {
		cw_tmp.back() ^= cw_tmp[i];
	}
	if (cw_tmp.back() != 0)
		return false;
	else
		return true;
}

