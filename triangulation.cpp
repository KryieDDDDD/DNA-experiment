#include "triangulation.h"

//global parameters generated and shared by triangulation
bool use_DSDS = true;//default
int num_rec;//number of recovered VNs (size of L)
int num_inact;//number of inactive VNs (size of R)
VI	CN_pos_idx;//CN pos (stores) idx
VI	VN_idx_pos;//VN idx (is stored at) pos
clock_t	time_triang;//clock time for triangulation

VI	CN_deg_tri;//upper-bound of CN degree
VI	CN_sum_tri;//sum of degrees of neighbouring VNs

VI	par_tri;//parent node of each VN
VI	rk_tri;//union heuristic: by rank
VI	wt_tri;//the weight of a rooted tree
int num_comp;//#component (size >= 2)

int tri_pq;

int find_root(int id) {//path compression: O(#find * alpha(#find, n))
	return	id == par_tri[id] ? id : par_tri[id] = find_root(par_tri[id]);
}

VI dexor(int CN_id) {
	VI	res;
	for (int i = 0; i < CN_adj[CN_id].size(); ++i) {
		if (VN_idx_pos[CN_adj[CN_id][i]] == -1) {//does not move to L or R
			res.push_back(CN_adj[CN_id][i]);
		}
	}
	return	res;
}

//定义优先队列排序规则；
//template <typename T>
class CMP
{
public:
	bool operator()(int CN_idx1, int CN_idx2) {
		//if (copy_EES[CN_idx1].back() != copy_EES[CN_idx2].back()) {
		//	return copy_EES[CN_idx1].back() < copy_EES[CN_idx2].back();
		//}
		//else {
		//	return	CN_adj[CN_idx1].size() < CN_adj[CN_idx2].size();
		//}
		return	CN_adj[CN_idx1].size() < CN_adj[CN_idx2].size();//注意优先队列和sort的默认排序规则相反，所以sort是>,而这里是<!
	}
};

VI recover(VI& CN_ids) {
	VI union_CNs;
	
	priority_queue<int, vector<int>, CMP>pq;
	for (int i = 0; i < CN_ids.size(); ++i) {
		pq.push(CN_ids[i]);
	}
	
	while (!CN_ids.empty()) {
		int CN_id;
		if (tri_pq == 1) {
			CN_id = pq.top();
			pq.pop();
			CN_ids.erase(remove(CN_ids.begin(), CN_ids.end(), CN_id), CN_ids.end());
		}
		else if (tri_pq == 0) {
			CN_id = CN_ids.back();
			CN_ids.pop_back();
		}
		//int CN_id = CN_ids.back();
		//CN_ids.pop_back();
		CN_pos_idx.push_back(CN_id);
		if (CN_deg_tri[CN_id] == 0) {
			continue;
		}
		else {
			assert(CN_deg_tri[CN_id] == 1);
		}
		
		VI cur_adj = dexor(CN_id);
		assert(cur_adj.size() == 1);
		
		int v = cur_adj[0];
		swap(CN_pos_idx.back(), CN_pos_idx[num_rec]);//make top rows form lower triangular matrix
		VN_idx_pos[v] = num_rec++;
		
		for (int i = 0; i < VN_adj[v].size(); ++i) {
			int c = VN_adj[v][i];
			--CN_deg_tri[c];
			CN_sum_tri[c] -= VN_adj[v].size();
			if (CN_deg_tri[c] == 1) {
				CN_ids.push_back(c);
				if (tri_pq == 1) {
					pq.push(c);
				}
			}
			else if (CN_deg_tri[c] == 2 && use_DSDS) {
				union_CNs.push_back(c);
			}
		}
	}
	
	return	union_CNs;
}

void union_root(VI& CN_ids) {
	while (!CN_ids.empty()) {
		int CN_id = CN_ids.back();
		CN_ids.pop_back();
		if (CN_deg_tri[CN_id] == 0) {
			continue;
		}
		else {
			assert(CN_deg_tri[CN_id] == 2);
		}
		
		VI cur_adj = dexor(CN_id);
		assert(cur_adj.size() == 2);
		
		int a = find_root(cur_adj[0]);
		int b = find_root(cur_adj[1]);
		if (a == b) {
			continue;
		}
		if (rk_tri[a] < rk_tri[b]) {
			swap(a, b);
		}
		par_tri[b] = a;
		rk_tri[a] += rk_tri[a] == rk_tri[b];
		num_comp += (wt_tri[a] == 1) + (wt_tri[b] == 1) - 1;
		wt_tri[a] += wt_tri[b];
	}
}

VI inactivation() {
	VI res;
	if (num_comp > 0) {
		for (int j = 0; j < par_tri.size(); ++j) {
			if (par_tri[j] != j || wt_tri[j] <= 1 || VN_idx_pos[j] != -1) {
				continue;
			}
			else if (res.empty()) {
				res.push_back(j);
			}
			else if (wt_tri[j] > wt_tri[res[0]]) {
				res[0] = j;
			}
		}
	}
	
	if (res.empty()) {//no degree-two row if use_DSDS
		int CN_id = -1;
		for (int i = 0; i < CN_deg_tri.size(); ++i) {
			if (CN_deg_tri[i] <= 1) {
				continue;
			}
			else if (CN_id == -1) {
				CN_id = i;
			}
			else if (CN_deg_tri[i] < CN_deg_tri[CN_id]) {
				CN_id = i;
			}
			else if (CN_deg_tri[i] == CN_deg_tri[CN_id] && CN_sum_tri[i] > CN_sum_tri[CN_id]) {
				CN_id = i;
			}
		}
		
		if (CN_id != -1) {
			res = dexor(CN_id);
			assert(CN_deg_tri[CN_id] == res.size());
			res.erase(res.begin());
		}
	}
	return	res;
}

//input: m, n, CN_adj
//get form [L R], where L is a lower triangular matrix, |L| = num_rec, |R| = num_inact
//permutation on rows/columns are recorded by CN_pos_idx and VN_idx_pos
//return true: num_rec + num_inact = n; false: otherwise
bool triangulation() {
	time_triang = clock();
	
	assert(m == CN_adj.size());
	VN_adj.clear();
	VN_adj.resize(n);
	CN_deg_tri.resize(m);
	
	VI CN1, CN2;
	for (int i = 0; i < m; ++i) {
		CN_deg_tri[i] = CN_adj[i].size();
		for (int j = 0; j < CN_adj[i].size(); ++j) {
			assert(CN_adj[i][j] >= 0 && CN_adj[i][j] < n);
			VN_adj[CN_adj[i][j]].push_back(i);
		}
		
		//allow CN_deg_tri[i] == 0 for LDPC code
		if (CN_deg_tri[i] == 1) {//a new ripple
			CN1.push_back(i);
		}
		else if (CN_deg_tri[i] == 2 && use_DSDS) {//a CN connecting two VNs
			CN2.push_back(i);
		}
	}
	
	CN_sum_tri.resize(m);
	for (int i = 0; i < m; ++i) {
		CN_sum_tri[i] = 0;
		for (int j = 0; j < CN_adj[i].size(); ++j) {
			CN_sum_tri[i] += VN_adj[CN_adj[i][j]].size();
		}
	}
	
	if (use_DSDS) {
		par_tri.resize(n);
		for (int i = 0; i < par_tri.size(); ++i) {
			par_tri[i] = i;
		}
		rk_tri.resize(n);
		fill(rk_tri.begin(), rk_tri.end(), 1);
		wt_tri.resize(n);
		fill(wt_tri.begin(), wt_tri.end(), 1);
	}
	else {
		par_tri.clear();//inactivation does not consider component
	}
	
	num_comp = 0;
	num_inact = 0;
	num_rec = 0;
	CN_pos_idx.clear();
	VN_idx_pos.resize(n);
	fill(VN_idx_pos.begin(), VN_idx_pos.end(), -1);
	
	while (true) {
		VI union_CNs = recover(CN1);//CN1 becomes empty
		union_root(union_CNs);
		union_root(CN2);//CN2 becomes empty
		//cout << num_rec << " + " << num_inact << " =? " << n << endl;
		if (num_rec + num_inact == n) {
			break;
		}
		
		VI inVNs = inactivation();
		if (inVNs.empty()) {
			time_triang = clock() - time_triang;
			return	false;//such as existing zero columns 
		}
		
		for (int j = 0; j < inVNs.size(); ++j) {
			int v = inVNs[j];
			++num_inact;
			VN_idx_pos[v] = n - num_inact;
			
			for (int i = 0; i < VN_adj[v].size(); ++i) {
				int c = VN_adj[v][i];
				CN_sum_tri[c] -= VN_adj[v].size();
				--CN_deg_tri[c];
				if (CN_deg_tri[c] == 1) {
					CN1.push_back(c);
				}
				else if (CN_deg_tri[c] == 2 && use_DSDS) {
					CN2.push_back(c);
				}
			}
		}
		
	}
	
	time_triang = clock() - time_triang;
	//at this point, first num_rec rows and columns are lower triangular matrix
	return	true;
}

