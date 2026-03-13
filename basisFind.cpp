#include "basisFind.h"

VB	CN_err;//1: CN is incorrect; 0: CN is correct
//int	times_n_cor_basis = 0;
//LL	num_inc_basis = 0;
//double	basis_wt = 0.0;
//VD	p_att_mean(2, 0.0);
//VD	p_att_var(2, 0.0);

//double CON_basis_P = 0.0;
//int CON_basis_N;

//input: m, n, l, CN_adj, CN_val, num_inact, CN_pos_idx, VN_idx_pos
//return 0: n correct basis elements and they uniquely attend top n linear representations; 
//1: error: less than n correct basis elements;
//2: error: n correct basis elements but they do not uniquely attend top n linear representations
int BFA_straightforward() {
	VI	num_basis(2, 0);//0/1: number of correct/incorrect basis elements
	VI	num_LR;//#LR (linear representation) each basis element attends
	VB	basis_err;
	VVI	B(n + l), Q(n + l);
	//int	sum_wt = 0;
	
	//CON_basis_N = 0;
	
	for (int i = 0;; ++i) {
		if (num_basis[1] > l || num_basis[0] + m - i < n) {
			return	1;//error due to correct basis is not full rank
		}
		else if (i == m) {
			break;
		}
		
		int min_v = -1;
		VI b = CN_val[CN_pos_idx[i]];
		b.resize((l + num_inact + BWI - 1) / BWI, 0);
		VI q((num_LR.size() + BWI - 1) / BWI, 0);
		for (int j = 0; j < CN_adj[CN_pos_idx[i]].size(); ++j) {
			assert(CN_adj[CN_pos_idx[i]][j] >= 0 && CN_adj[CN_pos_idx[i]][j] < n);
			int v = VN_idx_pos[CN_adj[CN_pos_idx[i]][j]];
			if (v + num_inact < n) {//v belongs to L
				if (B[v].empty()) {
					//cout << "v��ֵ��" << v << endl;
					//cout << "min_v��ֵ��" << min_v << endl;
					//cout << VN_idx_pos[CN_adj[CN_pos_idx[0]][0]] << "\t" << VN_idx_pos[CN_adj[CN_pos_idx[0]][1]] << endl;
					assert(min_v == -1);
					min_v = v;
					//break;
				}
				else {
					xorEq(b, B[v]);
					xorEq(q, Q[v]);//q may be longer than Q[v] 
				}
			}
			else {
				b[(n + l - 1 - v) / BWI] ^= (1 << (n + l - 1 - v) % BWI);
			}
		}
		
		for (int j = b.size() - 1; j >= 0 && min_v == -1; --j) {
			for (int t = BWI - 1; t >= 0 && b[j] != 0; --t) {
				if (b[j] & (1 << t)) {
					int cur_v = n + l - 1 - j * BWI - t;
					assert(cur_v + num_inact >= n);//cur_v belongs to R
					if (B[cur_v].empty()) {
						min_v = cur_v;
						break;
					}
					else {
						xorEq(b, B[cur_v]);
						xorEq(q, Q[cur_v]);
					}
				}
			}
		}
		
		if (min_v != -1) {//add into basis
			B[min_v] = b;
			Q[min_v] = q;
			if (num_LR.size() % BWI == 0) {
				Q[min_v].push_back(0);
			}
			Q[min_v].back() ^= (1 << num_LR.size() % BWI);
			//num_LR.push_back(0);
			
			
			//if (CN_con[CN_pos_idx[i]] == 1) {
			//	num_LR.push_back(m + 1);
			//	//CON_basis_N++;
			//}
			//else {
			//	num_LR.push_back(0);
			//}
			
			num_LR.push_back(copy_EES[CN_pos_idx[i]].back());
			basis_err.push_back(CN_err[CN_pos_idx[i]]);
			++num_basis[basis_err.back()];
			//sum_wt += CN_adj[CN_pos_idx[i]].size();
		}
		else {//b is a LC of the current basis; q indicates the basis elements
			for (int t = 0; t < q.size(); ++t) {
				if (q[t] != 0) {
					for (int j = 0; j < BWI; ++j) {
						if ((q[t] >> j) & 1) {
							num_LR[t * BWI + j] += copy_EES[CN_pos_idx[i]].back();
						}
					}
				}
			}
		}
	}
	//cout << num_basis[0] << "\t" << num_basis[1] << endl;
	//assert(1 == 2);
	
	assert(num_basis[0] == n && num_basis[1] <= l);//Ϊʲônum_basis[0]������>n����ΪB[v].size()�ܵ�VN_idx_pos.size()�����ơ�
	//++times_n_cor_basis;
	//num_inc_basis += num_basis[1];
	//basis_wt += 1.0 * sum_wt / (num_basis[0] + num_basis[1]);
	
	//ͳ��Լ�����ڻ��������ƽ�������Լ�����Լ���ж��ڻ�����ĸ��ʣ�
	//	CON_basis_P += 1.0 * CON_basis_N / (num_basis[0] + num_basis[1]);
	//	if (CON_basis_N == S + H) {
	//		CON_basis_AP++;
	//	}
	
	int cor_min = m;//minimum LRs of correct basis
	int err_max = -1;//maximum LRs of erroneous basis
	//	VLL	sum_LR(2, 0);
	//	VLL	sum2_LR(2, 0);
	for (int i = 0; i < num_LR.size(); ++i) {
		if (basis_err[i]) {
			err_max = max(err_max, num_LR[i]);
		}
		else {
			cor_min = min(cor_min, num_LR[i]);
		}
		
		//		sum_LR[basis_err[i]] += num_LR[i];
		//		sum2_LR[basis_err[i]] += 1LL * num_LR[i] * num_LR[i];
		
	}
	//if(cor_min < err_max || cor_min==err_max){
	//	for (int i = 0; i <num_LR.size(); ++i) {
	//		if (basis_err[i]) {
	//			cout << "��" << i << "����(����)�μ�LR������" << num_LR[i] << endl;
	//		}
	//		else {
	//			cout << "��" << i << "����(��ȷ)�μ�LR������" << num_LR[i] << endl;
	//		}
	//	}
	//	cout <<"������������LR������" <<err_max << "\t" <<"��ȷ��������СLR������" <<cor_min << endl;
	//	//cout << 1 << endl;
	//	//assert(1 == 2);
	//	cout << endl;
	//}
	
	//	double redunSym = max(1, m - num_basis[0] - num_basis[1]);
	//	for (int i = 0; i < 2; ++i) {
	//		double ave = sum_LR[i] / redunSym / max(num_basis[i], 1);
	//		p_att_mean[i] += ave;
	//		p_att_var[i] += sum2_LR[i] / redunSym / redunSym / max(num_basis[i], 1) - ave * ave;
	//	}
	
	return	cor_min > err_max ? 0 : 2;
}

bool cmp(int CN_idx1, int CN_idx2) {
	//return	CN_adj[CN_idx1].size() > CN_adj[CN_idx2].size();
	return pq_key[CN_idx1] > pq_key[CN_idx2];
}

//bool cmp1(int CN_idx1, int CN_idx2){
//	return abs(read_EES[CN_idx1].back()) > abs(read_EES[CN_idx2].back());
//}

//input: m, n, l, CN_adj, CN_val
//return 0: n correct basis elements and they uniquely attend top n linear representations; 
//1: error: less than n correct basis elements;
//2: error: n correct basis elements but they do not uniquely attend top n linear representations
int	BFA() {
	assert(m == CN_adj.size() && m == CN_val.size() && n > 0 && (l + BWI - 1) / BWI == CN_val[0].size());
	if (dec_type == "tri") {//triangulation-based BFA
		//get form [L R], where L is a lower triangular matrix; |R| = num_inact
		tri_pq = 1;
		if (!triangulation()) {
			return	1;//error due to existing zero columns
		}
		else {
			int s = num_rec;
			//num_inact = n;
			VI CN_pos_idx_1 = CN_pos_idx;
			//VI VN_idx_pos_1 = VN_idx_pos;
			CN_pos_idx.clear();
			CN_pos_idx.resize(m);
			//CN_pos_idx.resize(m);
			//VN_idx_pos.clear();
			//VN_idx_pos.resize(n);
			//cout << CN_pos_idx_1.size() << endl;
			//cout << endl;
			for (int i = 0; i < num_rec; ++i) {
				CN_pos_idx[i] = CN_pos_idx_1[i];
			}
			for (int i = num_rec; i < m; ++i) {
				if (CN_con[CN_pos_idx_1[i]] == 1) {
					CN_pos_idx[s] = CN_pos_idx_1[i];
					s++;
					//VN_idx_pos[VN_pos_idx[CN_idx_pos[CN_pos_idx_1[i]]]] = CN_idx_pos[CN_pos_idx_1[i]];
					//cout << "CN" << i << ": " << CN_pos_idx_1[i] << "\t" << "VN" << VN_pos_idx[CN_idx_pos[CN_pos_idx_1[i]]] << ": " << VN_idx_pos[VN_pos_idx[CN_idx_pos[CN_pos_idx_1[i]]]] << endl;
					//if (CN_pair_VN[CN_pos_idx_1[i]] == 1) {
					//	
					//}
				}
			}
			for (int i = num_rec; i < m; ++i) {
				if (CN_con[CN_pos_idx_1[i]] == 0) {
					CN_pos_idx[s] = CN_pos_idx_1[i];
					s++;
				}
			}
			//cout << endl;
			//cout << CN_pos_idx.size() << endl;
			//for (int i = 0; i < m; ++i) {
			//	if (!be_zeros(CN_val[CN_pos_idx_1[i]])) {
			//		CN_pos_idx.push_back(CN_pos_idx_1[i]);
			//		//VN_idx_pos[VN_pos_idx[CN_idx_pos[CN_pos_idx_1[i]]]] = CN_idx_pos[CN_pos_idx_1[i]];
			//		//cout << "CN" << i << ": " << CN_pos_idx_1[i] << "\t" << "VN" << VN_pos_idx[CN_idx_pos[CN_pos_idx_1[i]]] << ": " << VN_idx_pos[VN_pos_idx[CN_idx_pos[CN_pos_idx_1[i]]]] << endl;
			//	}
			//	//	if (CN_pair_VN[CN_pos_idx_1[i]] == 1) {
			//	//		
			//	//}
			//}
			/*VN_idx_pos.resize(n);
			for (int i = 0; i < VN_idx_pos.size(); ++i) {
				VN_idx_pos[i] = i;
			}
			assert(1 == 2);
			cout << endl;
			cout << VN_idx_pos.size() << endl;

			cout << CN_pos_idx[S+H-1] << endl;
			assert(CN_pos_idx.size() == S + H);
			for (int i = 0; i < m - S - H; ++i) {
				CN_pos_idx.push_back(CN_pos_idx_1[i]);
			}*/
			
		}
	}
	else {
		num_inact = n;//indicate L is empty and R is full
		if (dec_type == "sort") {
			for(int i = S + H; i < m;){
				int j = i;
				while(++ j < m && copy_EES[CN_pos_idx[j]].back() == copy_EES[CN_pos_idx[i]].back());
				if(j - i > 1){
					stable_sort(CN_pos_idx.begin() + i, CN_pos_idx.begin() + j, cmp);
				}
				i = j;
			}
		}
		
		//		CN_pos_idx.resize(m);
		//		for (int i = 0; i < m; ++i) {
		//			CN_pos_idx[i] = i;
		//		}
		
		//if (dec_type == "sort") {
		//	sort(CN_pos_idx.begin() + S + H, CN_pos_idx.end(), cmp);//sort received symbols from most weighted to least weighted
		//}
		
		//		VN_idx_pos.resize(n);
		//		for (int i = 0; i < VN_idx_pos.size(); ++i) {
		//			VN_idx_pos[i] = i;
		//		}
	}
	return	BFA_straightforward();
}

