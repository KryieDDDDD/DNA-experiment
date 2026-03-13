#include "decBEC.h"

bool use_BP = false;//default SGE
VVI rec_info;//recovered info: empty for non-recovered

VI CN_deg_dec;

bool recover_BP(VI& CN_ids) {
	int num_rec_BP = 0;
	while (!CN_ids.empty()) {
		int CN_id = CN_ids.back();
		CN_ids.pop_back();
		if (CN_deg_dec[CN_id] == 0) {
			continue;
		}
		else {
			assert(CN_deg_dec[CN_id] == 1);
		}
		
		VI val = CN_val[CN_id];
		int v = -1;
		for (int i = 0; i < CN_adj[CN_id].size(); ++i) {
			if (rec_info[CN_adj[CN_id][i]].empty()) {
				assert(v == -1);
				v = CN_adj[CN_id][i];
			}
			else {
				xorEq(val, rec_info[CN_adj[CN_id][i]]);
			}
		}
		assert(v >= 0);
		rec_info[v] = val;
		++num_rec_BP;
		
		for (int i = 0; i < VN_adj[v].size(); ++i) {
			if (--CN_deg_dec[VN_adj[v][i]] == 1) {
				CN_ids.push_back(VN_adj[v][i]);
			}
		}
	}
	
	return	num_rec_BP == n;
}

bool BP() {
	//cout << "BP here ----------" << endl;
	assert(m == CN_adj.size());
	rec_info.clear();
	rec_info.resize(n);
	VN_adj.clear();
	VN_adj.resize(n);
	CN_deg_dec.resize(m);
	
	VI CN1;
	for (int i = 0; i < m; ++i) {
		CN_deg_dec[i] = CN_adj[i].size();
		for (int j = 0; j < CN_adj[i].size(); ++j) {
			assert(CN_adj[i][j] >= 0 && CN_adj[i][j] < n);
			VN_adj[CN_adj[i][j]].push_back(i);
		}
		if (CN_deg_dec[i] == 1) {//a new ripple
			CN1.push_back(i);
		}
	}
	return	recover_BP(CN1);
}

//mat * x = val: return x or empty matrix for failure
//mat corresponds to a bool matrix; each integer correponds to [BWI] coefficients
VVI solve_eqn_by_GE(int num_unknown, VVI mat, VVI val) {
	assert(num_unknown > 0 && num_unknown <= mat.size() && mat.size() == val.size());
	if (mat.size() < num_unknown) {
		return	VVI(0, VI(0));
	}
	
	VI idx(num_unknown, -1);//idx[i]: idx of the row of mat whose leftmost 1 is at i
	int idx_size = 0;
	for (int i = 0; i < mat.size() && idx_size < num_unknown; ++i) {
		for (int j = 0; j < mat[i].size(); ++j) {
			for (int b = 0; b < BWI && mat[i][j] != 0; ++b) {
				if (mat[i][j] & (1 << b)) {//leftmost 1
					int id = j * BWI + b;
					if (idx[id] < 0) {//does not exist previously
						idx[id] = i;
						++idx_size;
						goto out_loop_j;
					}
					else {//leftmost 1 is eliminated
						xorEq(mat[i], mat[idx[id]]);
						xorEq(val[i], val[idx[id]]);
					}
				}
			}
		}
		out_loop_j:;
	}
	
	VVI res(num_unknown);
	for (int i = num_unknown - 1; i >= 0; --i) {
		if (idx[i] < 0) {
			return	VVI(0, VI(0));
		}
		
		res[i] = val[idx[i]];
		for (int j = i + 1; j < num_unknown; ++j) {
			if (mat[idx[i]][j / BWI] & (1 << j % BWI)) {
				xorEq(res[i], res[j]);
			}
		}
	}
	return	res;
}

bool SGE() {
	if (!triangulation()) {
		return	false;
	}
	else {
		assert(num_rec + num_inact == n);
	}
	
	rec_info.clear();
	rec_info.resize(n);
	
	if (num_rec < n) {
		VVI mat(m, VI((num_inact + BWI - 1) / BWI, 0));
		VVI val(m);
		for (int i = 0; i < m; ++i) {
			int c = CN_pos_idx[i];
			val[i] = CN_val[c];
			for (int j = 0; j < CN_adj[c].size(); ++j) {
				int v = VN_idx_pos[CN_adj[c][j]];
				assert(v >= 0 && v < n);
				
				if (v < num_rec && v < i) {
					xorEq(mat[i], mat[v]);
					xorEq(val[i], val[v]);
				}
				else if (v >= num_rec) {
					mat[i][(v - num_rec) / BWI] ^= (1 << (v - num_rec) % BWI);
				}
				else {
					assert(v == i);
				}
			}
		}
		
		mat.erase(mat.begin(), mat.begin() + num_rec);
		val.erase(val.begin(), val.begin() + num_rec);
		VVI rec_val = solve_eqn_by_GE(num_inact, mat, val);
		if (rec_val.empty()) {
			return	false;
		}
		else {
			for (int j = 0; j < n; ++j) {
				if (VN_idx_pos[j] >= num_rec) {
					rec_info[j] = rec_val[VN_idx_pos[j] - num_rec];
				}
			}
		}
	}
	
	for (int i = 0; i < num_rec; ++i) {
		int c = CN_pos_idx[i];
		int v = -1;
		VI val = CN_val[c];
		for (int j = 0; j < CN_adj[c].size(); ++j) {
			if (rec_info[CN_adj[c][j]].empty()) {
				assert(v == -1);
				v = CN_adj[c][j];
			}
			else {
				xorEq(val, rec_info[CN_adj[c][j]]);
			}
		}
		assert(v >= 0);
		rec_info[v] = val;
	}
	return	true;
}

//input: m, n, CN_adj, CN_val
//return true/false: rec_info has n/<n recovered source symbols
bool decBEC() {
	if (m < n) {
		return	false;
	}
	return	use_BP ? BP() : SGE();
}
