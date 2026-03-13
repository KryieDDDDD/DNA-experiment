//该程序得到国家重点研发计划项目(2021YFF1200200)的支持
//
//
#include "simulate.h"
int MAX_FE = 2;//maximum frame errors to stop decoding 
const int MIN_TC = 10;//minimum testcases to stop decoding
const int MAX_TC = 10e8;//minimum testcases to stop decoding

VI CN_con;
VVI read_EES;
string dec_type;//BP, tri, str, sort <-> BP, triangulation-based BFA, straightforward BFA, sorted-weight BFA
VVI copy_EES;
VD pq_key;

class L_dis {
public:
	int dis;
	int ins;
	int del;
	int sub;
	L_dis();
};
L_dis::L_dis() {
	dis = 0;
	ins = 0;
	del = 0;
	sub = 0;
}
int min_positive(int a, int b, int c) {
	VI temp;
	if (a >= 0)
		temp.push_back(a);
	if (b >= 0)
		temp.push_back(b);
	if (c >= 0)
		temp.push_back(c);
	if (temp.empty())
		return -2;
	else if (temp.size() == 1)
		return temp.back();
	else {
		int m = temp[0];
		for (int i = 1; i < temp.size(); ++i)
			if (m > temp[i])
				m = temp[i];
		return m;
	}
}
L_dis cal_L_distance(string s1, string s2) {//s2 -> s1
	L_dis L;
	int maxDis = 1;
	VVI dp(s1.size() + 1, VI((s2.size() + 1), -1));
	for (int i = 0; i < dp.size(); ++i)
		dp[i][0] = i;
	for (int j = 0; j < dp[0].size(); ++j)
		dp[0][j] = j;
	for (int i = 1; i < dp.size(); ++i) {
		//for (int j = 1; j < dp[i].size(); ++j) {
		for (int j = (i - maxDis <= 1) ? 1 : i - maxDis; j < dp[i].size() && j <= i + maxDis; ++j) {
			if (s1[i - 1] == s2[j - 1])
				dp[i][j] = dp[i - 1][j - 1];
			else {
				//int itemp = min_positive(dp[i][j - 1], dp[i - 1][j - 1], dp[i - 1][j]);
				//if (itemp == -1) {
				//	dp[i][j] = -1;
				//	cout << "-1.........................." << endl;
				//}
				//else
				//	dp[i][j] = itemp + 1;
				dp[i][j] = min_positive(dp[i][j - 1], dp[i - 1][j - 1], dp[i - 1][j]) + 1;
			}
		}
	}
	L.dis = dp.back().back();
	//if (L.dis < 0)
	//	return L;
	//if (L.dis > maxDis) {
	//	L.dis = -1;
	//	return L;
	//}

	return L;
}
bool cmp1(int CN_idx1, int CN_idx2) {
	return copy_EES[CN_idx1].back() > copy_EES[CN_idx2].back();
}

void set_distri_LT() {

	distri_LT.clear();
	distri_LT.resize(41, 0);
	distri_LT[1] = 0.00971;
	distri_LT[2] = 0.4580;
	distri_LT[3] = 0.2100;
	distri_LT[4] = 0.1130;
	distri_LT[10] = 0.1110;
	distri_LT[11] = 0.0797;
	distri_LT[40] = 0.0156;

	for (int i = 2; i < distri_LT.size(); ++i) {
		distri_LT[i] += distri_LT[i - 1];
	}
	//distri_LT[0] = 0;
	//double ave_wt = 0.0;
	for (int i = 1; i < distri_LT.size(); ++i) {
		distri_LT[i] /= distri_LT.back();
		//ave_wt += (distri_LT[i] - distri_LT[i - 1]) * i;
	}
	distri_LT[0] = -1;//never select 0
	distri_LT.back() += eps;//ensure distri_LT.back() > 1
	//cerr << "average weight of distri_LT: " << ave_wt << endl;
}

class R10_seed {
public:
	VI seed;
	int num;
	bool wrong;
	R10_seed();
};
R10_seed::R10_seed() {
	seed = VI(0);
	num = 0;
	wrong = 0;
}

struct VectorHasher {
	int operator()(const vector<int>& V) const {
		int hash = V.size();
		for (auto& i : V) {
			hash ^= i + 0x9e3779b9 + (hash << 6) + (hash >> 2);
		}
		return hash;
	}
};

vector<bool> a_from_b(int a, int b) {
	//cout << "sample num = " << a << endl;
	vector<bool> rst(b);
	for (int i = 0; i < a; ++i)
		rst[i] = 1;
	//srand(time(0));
	random_shuffle(rst.begin(), rst.end());
	return rst;
}

bool cmp_num_des(R10_seed sd1, R10_seed sd2) {
	return sd1.num > sd2.num;
}

bool cmp_num_des_vec(pair<VI, pair<bool, int>> p1, pair<VI, pair<bool, int>> p2) {
	return p1.second.second > p2.second.second;
}

bool cmp_num_des_vec_(pair<VI, int> p1, pair<VI, int> p2) {
	return abs(p1.second) > abs(p2.second);
	//return p1.second > p2.second;
}

void fountain_dec_from_txt() {
	//double allCoverage[] = { 4.2, 4.3, 4.4};
	n = 28500;
	double allCoverage[] = { 4.0, 4.1, 4.2, 4.3};
	int upper_bound = 0;
	int detect_only = 0;
	bool use_bfa = 1;/////////////////////////////////
	VS dirNames;
	for (int i = 0; i <= 3; ++i) {
		//if (i == 5)
		//	MAX_FE = 10;
		double avg = allCoverage[i];
		char dirName[100] = "";
		//sprintf(dirName, "scheme%d,BFATri,notDetect", scheme_no);/////////////////////////
		sprintf(dirName, "scheme1,results");/////////////////////////
		dirNames.push_back(dirName);
		string logFile;
		mkdir_(dirName);
		logFile = creatNewFile((dirNames.back() + "/log").c_str());
		//unsigned seed = getFileNum(logFile.c_str()) + (time(NULL) << 10);
		//unsigned seed = 0;

		ifstream ifs, ifs2;

		//ifs2.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\1_dec_wrong_idx.txt", ios::in);
		//if (!ifs.is_open() || !ifs2.is_open())
			//assert(0);

		int allLine = 11702066;
		int maxLine;
		if (detect_only)
			maxLine = 9836516;
		else
			maxLine = 11202021;
		int lineNum = 0;
		dec_type = "tri";
		R10_seed sd;
		vector<R10_seed> SD;
		unordered_map<VI, pair<bool, int>, VectorHasher> symbol_info;
		string stemp;
		int itemp;
		VI vitemp;
		int testCases = 0, err = 0, err_rank = 0;

		PRNG.seed(1);
		VVI info = rand_VVI(K, (l + BWI - 1) / BWI);
		inte.clear();
		if (l % BWI != 0) {
			int mask = (1 << l % BWI) - 1;
			for (int i = 0; i < info.size(); ++i) {
				info[i].back() &= mask;//make the tail bits be 0
			}
		}
		R10_encode(info);
		unsigned seed = getFileNum(logFile.c_str()) + \
			(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch()).count() << 10);
		//PRNG.seed(seed);
		srand(seed);
		cout << "seed = " << seed << endl;
		time_t start = clock();
		while ((err < MAX_FE || testCases < MIN_TC) && testCases < MAX_TC) {
			cout << "*************************" << endl;
			cout << "testCases = " << ++testCases << endl;
			//cout << "scheme : 0, SGE" << ", avg read = " << avg << endl;///////////////////
			printf("scheme = 1, upper_bound = %d, use_bfa = %d, detect_only = %d, coverage = %lf \n", use_bfa, upper_bound, detect_only, avg);
			printf("max run = %d, n = %d \n", max_run, n - S - H);

			lineNum = 0;
			time_t startTime = clock();
			VB sampling = a_from_b(30000 * avg, allLine);
			sampling.resize(maxLine);
			SD.clear();
			symbol_info.clear();
			//ifs.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\0_dec_with_idx.txt", ios::in);//////////////////////
			if (!detect_only)
				ifs.open("Experiment_Data\\scheme1\\s1_inner_decoding_results.txt", ios::in);//////////////////////
			else
				ifs.open("Experiment_Data\\scheme1\\s1_inner_detection_results.txt", ios::in);//////////////////////////
			if (!ifs.is_open())
				assert(0);
			while (getline(ifs, stemp)) {
				if (!sampling[lineNum++])
					continue;
				stringstream ss;
				ss << stemp;
				bool wrong;
				ss >> wrong;
				vitemp.clear();
				for (int i = 0; i < 16; ++i) {
					ss >> itemp;
					vitemp.push_back(itemp);
				}
				auto iter = symbol_info.find(vitemp);
				if (iter == symbol_info.end())
					symbol_info.insert({ vitemp, { wrong, 1 } });
				else {
					++iter->second.second;
				}
			}
			ifs.close();

			//unordered_map<int, int>temp;
			//for (auto s : symbol_info) {
			//	auto iter = temp.find(s.first[0]);
			//	if (iter != temp.end())
			//		temp[iter->first] = 2;
			//	else {
			//		if (s.second.first && s.second.second == 1)
			//			temp[s.first[0]] = -1;
			//		else if (!s.second.first && s.second.second == 1)
			//			temp[s.first[0]] = 1;
			//		else
			//			temp[s.first[0]] = 2;
			//	}
			//}
			//int numOf2 = 0;
			//for (auto t : temp)
			//	if (t.second == 2)
			//		++numOf2;
			//cout << "numOf2 = " << numOf2 << endl;
			//int unCoveredNumbers = 0, oneCopy = 0;
			//for (auto t : temp) {
			//	if (t.second == -1)
			//		++unCoveredNumbers;
			//	if (abs(t.second) == 1)
			//		++oneCopy;
			//}
			//cout << "unC1 = " << unCoveredNumbers << ", unC0 = " << (30000 - temp.size()) << ", one copy num = " << oneCopy << endl;
			//continue;


			unordered_map<int, pair<VI, pair<bool, int>>> symbol_info_unique;

			unordered_map<int, int> maxSameCopy;
			for (const auto& i : symbol_info) {
				auto iter = symbol_info_unique.find(i.first[0]);
				if (iter == symbol_info_unique.end())
					symbol_info_unique[i.first[0]] = i;
				else {

					if (i.second.second > iter->second.second.second)
						iter->second = i;
					else if (i.second.second == iter->second.second.second) {

						maxSameCopy[i.first[0]] = iter->second.second.second;
					}
				}
			}
			for (auto m : maxSameCopy) {
				auto iter = symbol_info_unique.find(m.first);
				if (iter == symbol_info_unique.end()) {
					cerr << "something wrong here !!" << endl;
					return;
				}
				if (iter->second.second.second <= m.second)
					symbol_info_unique.erase(m.first);
			}
			maxSameCopy.clear();


			vector<pair<VI, pair<bool, int>>> tmp_info;
			for (const auto& i : symbol_info_unique)///////////////////////////////////
				tmp_info.push_back(i.second);/////////////////////////////////////////////
			symbol_info.clear();
			symbol_info_unique.clear();
			cout << "After cluster, SD.size() == " << tmp_info.size() << endl;







			sort(tmp_info.begin(), tmp_info.end(), cmp_num_des_vec);
			int maxFreq = tmp_info[0].second.second;

			CN_adj.clear();
			CN_val.clear();

			CN_con.clear();
			CN_err.clear();

			for (int i = 0; i < S + H; ++i) {
				CN_adj.push_back(R10_LH_adj[i]);
				CN_val.push_back(R10_LH_val[i]);
				CN_con.push_back(1);
				CN_err.push_back(0);
			}
			if (!use_bfa) {
				int sz = 0;
				while (sz < tmp_info.size() && !tmp_info[sz].second.first)
					++sz;
				tmp_info.resize(sz);

			}


			for (int i = 0; i < tmp_info.size(); ++i) {
				CN_adj.push_back(R10_xor_VN(tmp_info[i].first[0]));
				CN_val.push_back(VI(tmp_info[i].first.begin() + 1, tmp_info[i].first.end()));
				CN_con.push_back(0);
				CN_err.push_back(tmp_info[i].second.first);
			}
			cout << "CN_adj.size() == " << CN_adj.size() << endl;

			n = K + S + H;
			m = tmp_info.size() + S + H;

			if (!use_bfa) {
				cout << "use_bp = " << use_BP << ", start decBEC ++++++++++++++" << endl;
				if (!decBEC())
					++err;
			}
			else {//perform BFA to decode raptor code
				startTime = clock();
				CN_pos_idx.resize(m);
				for (int i = 0; i < m; ++i) {
					CN_pos_idx[i] = i;
				}
				VN_idx_pos.resize(n);
				for (int i = 0; i < VN_idx_pos.size(); ++i) {
					VN_idx_pos[i] = i;
				}

				copy_EES.clear();
				copy_EES.resize(m);
				for (int i = 0; i < S + H; ++i) {
					for (int j = 0; j < tmp_info[0].first.size(); ++j) {
						copy_EES[i].push_back(1);
					}
					copy_EES[i].push_back(m + 1);
				}
				for (int i = S + H; i < m; ++i) {
					for (int j = 0; j < tmp_info[0].first.size(); ++j) {
						copy_EES[i].push_back(tmp_info[i - S - H].first[j]);
					}
					//copy_EES[i].push_back(tmp_info[i - S - H].second.second);
					if (upper_bound)
						copy_EES[i].push_back(tmp_info[i - S - H].second.second);
					else 
						copy_EES[i].push_back(abs(tmp_info[i - S - H].second.second));
					//if (!upper_bound) {
					//	if (tmp_info[i - S - H].second.first)
					//		copy_EES[i].push_back(-tmp_info[i - S - H].second.second);
					//	else
					//		copy_EES[i].push_back(tmp_info[i - S - H].second.second);
					//}
				}

				stable_sort(CN_pos_idx.begin(), CN_pos_idx.end(), cmp1);

				cout << "start BFA --------------" << endl;
				//if (!decBEC()) {
				//	cout << "SGE unsuccess...." << endl;
				//	++err;
				//}
				//else
				//	cout << "SGE success!!!!" << endl;

				pq_key.resize(m);
				for (int i = 0; i < S + H; ++i)
					pq_key[i] = 64;
				for (int i = S + H; i < m; i++) {
					pq_key[i] = 32.0 * tmp_info[i - S - H].second.second / maxFreq + 1.0 * CN_adj[i].size() / n;
				}

				int res = BFA();
				if (res == 1) {
					++err;
					++err_rank;
				}
				else if (res == 2) {
					++err;
				}
				cout << "time 4 BFA = " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
			}
			cout << "err : " << err << "(" << err_rank << ")" << " of " << testCases << ", " << 1.0 * err / testCases << endl;
			if (clock() > start + 10 * CLOCKS_PER_SEC || testCases == MIN_TC \
				|| testCases == MAX_TC || err == MAX_FE) {
				start = clock();
				FILE* fp = fopen(logFile.c_str(), "w");
				fprintf(fp, "seed              %u\n\n", seed);
				fprintf(fp, "coverage          %lf\n", avg);
				fprintf(fp, "testCases         %d\n", testCases);
				fprintf(fp, "wordErr           %d\n", err);
				fprintf(fp, "FER               %f\n", 1.0 * err / testCases);
				fclose(fp);
			}
		}

	}
}


void SPCC_N_foutain_dec_from_txt() {
	//ifstream ifs1, ifs2;
	//cout << "1" << endl;
	n = 28970;
	ofstream ofs1, ofs2;
	//ifs1.open("C:\\Users\\xinmatrix\\Desktop\\Data\\1.txt", ios::in);
	//ifs2.open("C:\\Users\\xinmatrix\\Desktop\\Data\\C1.txt", ios::in);
	double allCoverage[] = { 5.6, 5.7, 5.8, 5.9, 15, 12 };
	//double allCoverage[] = {4, 4.5, 5, 5.5, 6};
	int scheme_no = 2;
	max_run = 3;
	bool use_bfa = 1;
	bool detect_only = 1;
	VS dirNames;
	ifstream ifs2;
	ifs2.open("DNA_experiment\\Experiment_Data\\scheme2\\s2_R10_symbols.txt", ios::in);////////////////////
	if (!ifs2.is_open()) {
		cout << "file not open!" << endl;
		assert(0);
	}
	string stemp;
	int itemp, itemp2;
	char ctemp;
	VI vitemp;
	unordered_map<int, VI> original_strand;
	while (getline(ifs2, stemp)) {
		if (stemp.size() == 1)
			continue;
		else if (stemp[stemp.size() - 1] == '\r')
			stemp.resize(stemp.size() - 1);
		stringstream ss;
		ss << stemp;
		ss >> itemp;
		vitemp.clear();
		vitemp.push_back(itemp);
		for (int i = 0; i < 15; ++i) {
			ss >> itemp2;
			vitemp.push_back(itemp2);
		}
		original_strand[itemp] = vitemp;
	}

	//ofstream ofss;
	//ofss.open("/mnt/dingyi/testInnerCode.txt", ios::out);
	int maxErr = 100;
	for (int i = 0; i <= 3; ++i) {

		//if (i == 3)
		//	MAX_FE = 10;
		double avg = allCoverage[i];
		char dirName[100] = "";
		sprintf(dirName, "scheme%d,results", scheme_no);///////////////////////////////
		dirNames.push_back(dirName);
		string logFile;
		mkdir_(dirName);
		logFile = creatNewFile((dirNames.back() + "/log").c_str());
		unsigned seed = getFileNum(logFile.c_str()) + \
			(chrono::duration_cast<chrono::microseconds>(chrono::system_clock::now().time_since_epoch()).count() << 10);
		srand(seed);
		cout << "seed = " << seed << endl;
		ifstream ifs1;
		//ifs2.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\2_R10_Symbol.txt", ios::in);////////////////////
		//ofs1.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\2_dec.txt", ios::out);
		//ofs2.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\2_decInfo.txt", ios::out);

		//if (!ifs1.is_open() || !ifs2.is_open() || !ofs1.is_open() || !ofs2.is_open())
			//assert(0);

		int t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0;//t1 - no error; t2 - dec suc; t3 - undetec error; t4 - no dec res; t5 - < 249 || > 251
		long long dis = 0;

		//K = 28970;
		//n = 28970;

		int maxLine = 12538257;//2-12538257,3-11916498
		int maxLine_detect = 10597497;//2_detect_only.txt
		VVI record(50, VI(6));
		VI record_(3);

		unordered_map<string, unordered_map<string, int>> received_strand;
		if (K != 28970)
			cout << "K != 28970 !!!!!!" << endl;
		dec_type = "sort";

		//srand(0);
		unordered_map<VI, int, VectorHasher> R10_symbol;
		//string stemp;
		int itemp, itemp2;
		VI vitemp;
		int testCases = 0, err = 0, err_rank = 0;
		int aa = 0;
		int bb = 0;


		PRNG.seed(scheme_no);/////////////////////////////
		VVI info = rand_VVI(K, (l + BWI - 1) / BWI);
		inte.clear();
		if (l % BWI != 0) {
			int mask = (1 << l % BWI) - 1;
			for (int i = 0; i < info.size(); ++i) {
				info[i].back() &= mask;//make the tail bits be 0
			}
		}
		R10_encode(info);
		for (int i = 0; i < S + H; ++i) {
			CN_adj.push_back(R10_LH_adj[i]);
			CN_val.push_back(R10_LH_val[i]);
			CN_con.push_back(1);
			CN_err.push_back(0);
		}

		time_t start = clock();

		while ((err < MAX_FE || testCases < MIN_TC) && testCases < MAX_TC) {
			cout << "*************************" << endl;
			cout << "testCases = " << ++testCases << endl;
			//cout << "scheme : 2, SGE" << ", avg read = " << avg << endl;/////////////
			printf("scheme = %d, use_bfa = %d, detect_only = %d, coverage = %lf\n", scheme_no, use_bfa, detect_only, avg);
			//ifs1.open("C:\\Users\\xinmatrix\\Desktop\\Data\\2.txt", ios::in);////////////////////////
			//if (detect_only)
			//	ifs1.open("DNA_experiment\\Experiment_Data\\scheme2\\s2_detect_only.txt", ios::in);
			//else
			ifs1.open("DNA_experiment\\Experiment_Data\\scheme2\\s2_sequencing_oligos.txt", ios::in);////////////////////////
			if (!ifs1.is_open())
				assert(0);
			int lineNum = 0, decNum = 0;
			VB sampling = a_from_b(30000 * avg, maxLine);/////////////12538257, 11916498
			//if (detect_only)
			//	sampling.resize(maxLine_detect);
			received_strand.clear();
			time_t startTime = clock();
			itemp = 0;
			itemp2 = 0;
			while (getline(ifs1, stemp)) {
				if (!sampling[lineNum++])
					continue;
				if (stemp.size() == 1)
					continue;
				else if (stemp[stemp.size() - 1] == '\r')
					stemp.resize(stemp.size() - 1);
				 if (detect_only) {
				 	vitemp = bin2qua(string2bin(string(stemp.begin() + 1, stemp.end() - 1)));
				 	if (!valid_SPCC_cw(vitemp)) {
				 		continue;
				 	}
				 }
				string seed = string(stemp.begin(), stemp.begin() + 9);
				auto iter = received_strand.find(seed);
				if (iter == received_strand.end())
					received_strand[seed][stemp] = 1;
				else {

					auto iter2 = iter->second.find(stemp);
					if (iter2 == iter->second.end())
						iter->second[stemp] = 1;
					else
						//iter->second[stemp] = iter->second[stemp] + 1;
						++iter2->second;
				}

			}


			ifs1.close();
			cout << "after clustering, cluster size = " << received_strand.size() << ", time consuming : " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
			startTime = clock();
			lineNum = 0;
			//while (getline(ifs1, stemp)) {
			R10_symbol.clear();
			for (auto RS : received_strand) {
				++lineNum;
				//cout << "  lineNUm = " << lineNum << " ";
				// if (!(lineNum % 10000)) {
				// 	cout << 1.0 * lineNum / received_strand.size() << ", " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
				// }
				// if (record_[1] > maxErr) {////////////////////////////
				// 	int SPCCCases = 0;//////////////////////////////////////////////////////////
				// 	for (auto r_ : record_)
				// 		SPCCCases += r_;

				// 	ofss << "coverage = " << avg << ", SPCC testCases = " << SPCCCases << endl;
				// 	ofss << "right = " << record_[0] << ", " << 1.0 * record_[0] / SPCCCases << endl;
				// 	ofss << "wrong = " << record_[1] << ", " << 1.0 * record_[1] / SPCCCases << endl;
				// 	ofss << "erasure = " << record_[2] << ", " << 1.0 * record_[2] / SPCCCases << endl;
				// 	ofss << endl << endl;
				// 	break;
				// }
				int reads = 0;
				for (const auto& a : RS.second)
					reads += a.second;
				//cout << reads << " "<< endl;
				if (reads < record.size())
					++record[reads][0];
				char ctemp = RS.first[0];
				VVB strand_tmp;
				for (auto rs : RS.second) {
					for (int i = 0; i < rs.second; ++i) {
						stemp = string(rs.first.begin() + 1, rs.first.end() - 1);
						strand_tmp.push_back(string2bin(stemp));
					}
				}
				//cout << "0";
				int freq;
				VB dec_bin = SPCC_dec_4A_cluster(strand_tmp, freq);
				//cout << "1";
				if (dec_bin.empty()) {
					//++t1;
					++record_[2];
					if (reads < record.size())
						++record[reads][1];
					continue;
				}
				stemp = bin2string(dec_bin);
				stemp = ctemp + stemp + ctemp;
				if (!valid_GC_run(string2bin(stemp))) {
					//++t2;
					++record_[1];
					//cout << "record_[1] = " << record_[1] << endl;
					if (reads < record.size())
						++record[reads][2];
					continue;
				}
				//cout << "2";
				dec_bin = SPCC_extract_info(dec_bin);
				VB seed = VB(dec_bin.begin(), dec_bin.begin() + bw_seed);
				reverse(seed.begin(), seed.end());
				dec_bin = VB(dec_bin.begin() + bw_seed, dec_bin.end());
				dec_bin = SRT_decode(dec_bin);
				if (dec_bin.empty()) {
					//++t6;
					//cout << "t6++, " << lineNum << endl;
					++record_[1];
					//cout << "record_[1] = " << record_[1] << endl;
					if (reads < record.size())
						++record[reads][3];
					continue;
				}
				//cout << "3";
				dec_bin.insert(dec_bin.begin(), seed.begin(), seed.end());
				VI dec = compress_ES(dec_bin);
				bool existFlag = 0;
				auto _iter = original_strand.find(dec[0]);
				if (_iter != original_strand.end())
					if (_iter->second == dec)
						existFlag = 1;
				//if (freq > reads)
				//	cout << freq << " > " << reads << endl;
				if (existFlag) {
					//++t3;
					++record_[0];
					if (reads < record.size())
						++record[reads][4];
				}
				else {
					++record_[1];
					//cout << "err = " << record_[1] << endl;
					++t4;
					if (reads < record.size())
						++record[reads][5];
				}
				// if (record_[1] > maxErr) {
				// 	int SPCCCases = 0;//////////////////////////////////////////////////////////
				// 	for (auto r_ : record_)
				// 		SPCCCases += r_;

				// 	ofss << "coverage = " << avg << ", SPCC testCases = " << SPCCCases << endl;
				// 	ofss << "right = " << record_[0] << ", " << 1.0 * record_[0] / SPCCCases << endl;
				// 	ofss << "wrong = " << record_[1] << ", " << 1.0 * record_[1] / SPCCCases << endl;
				// 	ofss << "erasure = " << record_[2] << ", " << 1.0 * record_[2] / SPCCCases << endl;
				// 	ofss << endl << endl;
				// 	break;
				// }
				//continue;/////////////////////////////////////////////////
				//cout << "4";
				//if (freq == 1 && reads > 1) {
				//	++aa;
				//	bb += reads;
				//}
				auto R10_iter = R10_symbol.find(dec);
				if (R10_iter != R10_symbol.end()) {
					//cout << "exist !!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
					if (R10_iter->second > 0)
						R10_iter->second += freq;
					else
						R10_iter->second -= freq;
				}
				else {
					if (existFlag)
						R10_symbol[dec] = freq;
					else
						R10_symbol[dec] = -freq;
				}
				// if (!existFlag && freq >= 3) {
				// 	cout << "lineNum = " << lineNum << endl;
				// 	for (int i = 0; i < dec.size(); ++i)
				// 		cout << dec[i] << " ";
				// 	cout << endl;
				// }
			}
			// if (record_[1] > maxErr)
			// 	break;
			// continue;////////////////////////////////


			received_strand.clear();
			vector<pair<VI, int>>tmp_symbol;
			for (const auto& i : R10_symbol)
				tmp_symbol.push_back(i);
			R10_symbol.clear();
			cout << "After cluster, SD.size() == " << tmp_symbol.size() << endl;
			sort(tmp_symbol.begin(), tmp_symbol.end(), cmp_num_des_vec_);//////////////////////////////
			int maxFreq = tmp_symbol.front().second;
			//cout << "aa = " << aa << endl;
			//cout << "bb = " << bb << endl;

			CN_adj.resize(S + H);
			CN_val.resize(S + H);

			CN_con.resize(S + H);
			CN_err.resize(S + H);

			if (!use_bfa) {
				int itmp = 0;
				//while (itmp < tmp_symbol.size() && abs(tmp_symbol[itmp].second) > 1)
				//	++itmp;
				//cout << "reads > 1 : " << itmp << ", " << 1.0 * itmp / tmp_symbol.size() << endl;
				int sz = 0;
				while (sz < tmp_symbol.size() && tmp_symbol[sz].second > 0)
					++sz;
				cout << tmp_symbol[sz - 1].second << " -> " << tmp_symbol[sz].second << " ";
				//cout << itmp - 1 << ": " << tmp_symbol[itmp - 1].second << endl;
				//cout << itmp << ": " << tmp_symbol[itmp].second << endl;
				//cout << sz - 1 << ": " << tmp_symbol[sz - 1].second << endl;
				//cout << sz << ": " << tmp_symbol[sz].second << endl;
				tmp_symbol.resize(sz);

			}
			for (int i = 0; i < tmp_symbol.size(); ++i) {
				CN_adj.push_back(R10_xor_VN(tmp_symbol[i].first[0]));
				CN_val.push_back(VI(tmp_symbol[i].first.begin() + 1, tmp_symbol[i].first.end()));
				CN_con.push_back(0);
				if (tmp_symbol[i].second > 0)
					CN_err.push_back(0);
				else
					CN_err.push_back(1);
			}
			cout << "CN_adj.size() == " << CN_adj.size() << endl;
			n = K + S + H;
			m = tmp_symbol.size() + S + H;
			if (!use_bfa) {
				cout << "start SGE ++++++++++++++" << endl;
				if (!decBEC())
					++err;
			}
			else {
				startTime = clock();
				CN_pos_idx.resize(m);
				for (int i = 0; i < m; ++i) {
					CN_pos_idx[i] = i;
				}
				VN_idx_pos.resize(n);
				for (int i = 0; i < VN_idx_pos.size(); ++i) {
					VN_idx_pos[i] = i;
				}

				copy_EES.clear();
				copy_EES.resize(m);
				for (int i = 0; i < S + H; ++i) {
					for (int j = 0; j < tmp_symbol[0].first.size(); ++j) {
						copy_EES[i].push_back(1);
					}
					copy_EES[i].push_back(m + 1);
				}
				for (int i = S + H; i < m; ++i) {
					for (int j = 0; j < tmp_symbol[0].first.size(); ++j) {
						copy_EES[i].push_back(tmp_symbol[i - S - H].first[j]);
					}
					copy_EES[i].push_back(abs(tmp_symbol[i - S - H].second));
				}

				stable_sort(CN_pos_idx.begin(), CN_pos_idx.end(), cmp1);
				//int right = 0, wrong = 0, erasure = 0;
				//cout << "wrong strand copy nums : ";
				//for (auto const& a : tmp_symbol) {
				//	if (a.second > 0)
				//		++right;
				//	else {
				//		++wrong;
				//		cout << a.second << " ";
				//	}
				//}
				//cout << endl;
				//erasure = 30000 - right;
				//cout << "right wrong erasure : " << right << " " << wrong << " " << erasure << endl;


				cout << "start BFA --------------" << endl;
				//if (!decBEC()) {
				//	cout << "SGE unsuccess...." << endl;
				//	++err;
				//}
				//else
				//	cout << "SGE success!!!!" << endl;

				pq_key.resize(m);
				for (int i = 0; i < S + H; ++i)
					pq_key[i] = 64;
				for (int i = S + H; i < m; i++) {
					//pq_key[i] = CN_adj[i].size();
					//itemp = addr[i] & ((1 << addrBW) - 1);
					pq_key[i] = 32.0 * abs(tmp_symbol[i - S - H].second) / maxFreq + 1.0 * CN_adj[i].size() / n;
				}

				int res = BFA();
				if (res == 1) {
					++err;
					++err_rank;
				}
				else if (res == 2) {
					++err;
				}
				//else
				//	cout << "BFA return : " << res << endl;

				cout << "time 4 BFA = " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
			}
			//cout << "error strand = " << t4 << endl;
			cout << "err : " << err << "(" << err_rank << ")" << " of " << testCases << ", " << 1.0 * err / testCases << endl;
			if (!(testCases % 50)) {
				cout << "inner code dec info 4 diff copys : " << endl;
				for (int i = 0; i < record.size(); ++i) {
					if (!record[i][0])
						continue;
					cout << i << ", " << record[i][0] << ": ";
					for (int j = 1; j < record[i].size(); ++j)
						cout << record[i][j] << ", " << 1.0 * record[i][j] / record[i][0] << "   ";
					cout << endl;
				}
			}
			if (clock() > start + 10 * CLOCKS_PER_SEC || testCases == MIN_TC \
				|| testCases == MAX_TC || err == MAX_FE) {
				start = clock();
				FILE* fp = fopen(logFile.c_str(), "w");
				fprintf(fp, "seed              %u\n\n", seed);
				fprintf(fp, "coverage          %lf\n", avg);
				fprintf(fp, "testCases         %d\n", testCases);
				fprintf(fp, "wordErr           %d\n", err);
				fprintf(fp, "FER               %f\n", 1.0 * err / testCases);
				fclose(fp);
			}

		}
		//cout << "inner code dec info 4 diff copys : " << endl;
		// for (int i = 0; i < record.size(); ++i) {
		// 	if (!record[i][0])
		// 		continue;
		// 	cout << i << ", " << record[i][0] << ": ";
		// 	for (int j = 1; j < record[i].size(); ++j)
		// 		cout << record[i][j] << ", " << 1.0 * record[i][j] / record[i][0] << "   ";
		// 	cout << endl;
		// }
	}
}

void CopyDistributionAfterDecoding2() {
	int a = 0;
	ifstream ifs2;
	ofstream ofs, ofs2;
	ofs.open("C:\\Users\\xinmatrix\\Desktop\\2copyDistriAfterDecoding.txt", ios::out);
	ifs2.open("C:\\Users\\xinmatrix\\Desktop\\Data\\C2.txt", ios::in);

	//if (!ifs2.is_open() || !ifs1.is_open())
		//assert(0);
	string stemp;
	char ctemp;
	int lineNum = 0;
	int maxLine1 = 12538257;
	int maxLine2 = 11916498;
	unordered_map<string, int> original_strand;
	while (getline(ifs2, stemp)) {
		original_strand[string(stemp.begin(), stemp.begin() + 9)] = 0;
	}
	int t0 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0;
	time_t startTime = clock();
	for (double coverage = 5; coverage <= 6; coverage += 0.1) {
		cout << "coverage = " << coverage << endl;
		lineNum = 0;
		for (auto& o : original_strand)
			original_strand[o.first] = 0;
		VB sampling = a_from_b(coverage * 30000, maxLine1);
		ifstream ifs1;
		ifs1.open("C:\\Users\\xinmatrix\\Desktop\\Data\\2.txt", ios::in);
		if (!ifs1.is_open())
			cerr << "file not open !!" << endl;
		while (getline(ifs1, stemp)) {
			//if (stemp.size() < 249 || stemp.size() > 251) {
			//	++t0;
			//	continue;
			//}

			if (!(lineNum % 500000)) {
				cout << 1.0 * lineNum / maxLine1 << ", " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
				//cout << "-------- " << t0 << ", " << t1 << ", " << t2 << endl;
			}
			if (!sampling[lineNum++])
				continue;
			VVB vvbtemp = SPCC_dec_4A_seq(string2bin(string(stemp.begin() + 1, stemp.end() - 1)));
			if (!vvbtemp.size())
				continue;
			string sd = string(stemp.begin(), stemp.begin() + 9);
			auto iter = original_strand.find(sd);
			if (iter != original_strand.end())
				++iter->second;

		}
		unordered_map<int, int> copyInfo;
		for (const auto& o : original_strand) {
			auto iter = copyInfo.find(o.second);
			if (iter != copyInfo.end())
				++iter->second;
			else
				copyInfo[o.second] = 1;
		}
		vector<pair<int, int>> copyInfoVector;
		for (const auto& c : copyInfo)
			copyInfoVector.push_back(c);



		for (int i = 0; i < copyInfoVector.size(); ++i)
			for (int j = i + 1; j < copyInfoVector.size(); ++j)
				if (copyInfoVector[i].first > copyInfoVector[j].first)
					swap(copyInfoVector[i], copyInfoVector[j]);
		ofs << "coverage = " << coverage << ", copy distribution after decoding: " << endl;
		for (int i = 0; i < copyInfoVector.size(); ++i)
			ofs << copyInfoVector[i].first << " " << copyInfoVector[i].second << " " << 1.0 * copyInfoVector[i].second / 30000 << endl;
	}
}



int cal_hanmming_dis(string s1, string s2) {//binary
	//assert(s1.size() == s2.size());
	vector<bool> vb1 = string2bin(s1);
	vector<bool> vb2 = string2bin(s2);

	int sz = vb1.size() <= vb2.size() ? vb1.size() : vb2.size();
	int rst = 0;
	if (vb1.size() > vb2.size())
		rst += vb1.size() - vb2.size();
	else
		rst += vb2.size() - vb1.size();
	rst *= 2;
	for (int i = 0; i < sz; ++i)
		if (vb1[i] != vb2[i])
			++rst;
	return rst;

}
int cal_hanmming_dis(VB vb1, VB vb2) {
	int sz = vb1.size() <= vb2.size() ? vb1.size() : vb2.size();
	int rst = 0;
	if (vb1.size() > vb2.size())
		rst += vb1.size() - vb2.size();
	else
		rst += vb2.size() - vb1.size();
	rst *= 2;
	for (int i = 0; i < sz; ++i)
		if (vb1[i] != vb2[i])
			++rst;
	return rst;
}

void write_file() {
	ifstream ifs1, ifs2;
	ofstream ofs;
	ifs1.open("c:\\users\\xinmatrix\\desktop\\data_dec\\2_R10_Symbol.txt", ios::in);
	ifs2.open("c:\\users\\xinmatrix\\desktop\\data_dec\\2_dec.txt", ios::in);
	ofs.open("c:\\users\\xinmatrix\\desktop\\data_dec\\2_dec_with_idx.txt", ios::out);
	if (!ifs1.is_open() || !ifs2.is_open() || !ofs.is_open())
		assert(0);
	string stemp;
	int maxLine = 11654831;
	int lineNum = 0;
	int itemp = 0, itemp2 = 0, right = 0, wrong = 0;
	bool exist = 0;
	VI vitemp;
	unordered_map<int, VI> symbol;
	while (getline(ifs1, stemp)) {
		stringstream ss;
		ss << stemp;
		ss >> itemp;
		vitemp.clear();
		for (int i = 0; i < 15; ++i) {
			ss >> itemp2;
			vitemp.push_back(itemp2);
		}
		symbol[itemp] = vitemp;
	}
	clock_t startTime = clock();
	while (getline(ifs2, stemp)) {
		++lineNum;
		if (!(lineNum % 1000000))
			cout << 1.0 * lineNum / maxLine << ", " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
		exist = 0;
		stringstream ss;
		ss << stemp;
		ss >> itemp;
		auto iter = symbol.find(itemp);
		vitemp.clear();
		if (iter != symbol.end()) {
			for (int i = 0; i < 15; ++i) {
				ss >> itemp2;
				vitemp.push_back(itemp2);
			}
			if (vitemp == iter->second)
				exist = 1;
		}
		ofs << !exist << " " << stemp << endl;
		if (exist)
			++right;
		else
			++wrong;

	}

	cout << "right = " << right << endl;
	cout << "wrong = " << wrong << endl;

}

void testRealCoverage() {
	int scheme_no = 0;
	int detect_only = 0;
	bool use_bfa = 1;/////////////////////////////////
	int allLine = (scheme_no == 0) ? 11702066 : 12173158;///////////////////////11702066, 12173158
	int maxLine;
	if (detect_only)
		maxLine = (scheme_no == 0) ? 9836516 : 10230959;//////////////////detect:9836516, 10230959
	else
		maxLine = (scheme_no == 0) ? 11202021 : 11654831;///////////////////////11202021, 11654831
	unsigned SEED = time(NULL);
	srand(SEED);
	ifstream ifs;
	string sTemp;
	VI viTemp;
	int iTemp;
	int lineNum = 0;

	unordered_map<int, int> receivedSymbols;
	//fountain_dec_from_txt();
	for (double i = 3; i += 0.1; i <= 4.5) {
		cout << "coverage = " << i << " ";
		VB sampling = a_from_b(30000 * i, allLine);
		sampling.resize(maxLine);
		ifs.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\0_dec_with_idx.txt", ios::in);//////////////////////
		lineNum = 0;
		receivedSymbols.clear();
		while (getline(ifs, sTemp)) {
			//if (!(lineNum % 5000000))
			//	cout << 1.0 * lineNum / maxLine << ", " << 1.0 * (clock() - startTime) / CLK_TCK << endl;
			if (!sampling[lineNum++])
				continue;

			stringstream ss;
			ss << sTemp;
			bool wrongFlag;
			ss >> wrongFlag;
			viTemp.clear();
			for (int i = 0; i < 16; ++i) {
				ss >> iTemp;
				viTemp.push_back(iTemp);
			}
			auto iterCluster = receivedSymbols.find(viTemp[0]);
			if (iterCluster == receivedSymbols.end()) {
				if (wrongFlag) viTemp.push_back(-1);
				else viTemp.push_back(1);
				receivedSymbols[viTemp[0]] = wrongFlag;
			}
			//symbol_info[vitemp] = pair<bool, int>(wrong, 1);
			else {
				receivedSymbols[viTemp[0]] = 2;
			}
		}
		int uncovered = 30000 - receivedSymbols.size();
		cout << ", size here = " << uncovered << " ";
		for (auto a : receivedSymbols) {
			if (a.second == 1 || a.second == 0)
				++uncovered;
		}
		cout << ", final size = " << uncovered << endl;
		ifs.close();
	}
}

void generate_detect_only_file() {
	ifstream ifs;
	ofstream ofs;
	ifs.open("/mnt/dingyi/DNA_exp/Data/2.txt", ios::in);////////////////////////
	ofs.open("/mnt/dingyi/DNA_exp/Data/2_detect_only.txt", ios::out);
	if (!ifs.is_open())
		assert(0);
	int lineNum = 0;
	int outputNum = 0;
	int maxLine = 12538257;
	string stemp;
	clock_t startTime = clock();
	while (getline(ifs, stemp)) {
		++lineNum;
		if (stemp.size() == 1)
			continue;
		else if (stemp[stemp.size() - 1] == '\r')
			stemp.resize(stemp.size() - 1);
		if (!(lineNum % 100000)) {
			cout << 1.0 * lineNum / maxLine << ", " << 1.0 * (clock() - startTime) / CLOCKS_PER_SEC << endl;
			cout << "output num = " << outputNum << endl;
		}
		VI vitemp = bin2qua(string2bin(string(stemp.begin() + 1, stemp.end() - 1)));
		if (valid_SPCC_cw(vitemp)) {
			++outputNum;
			ofs << stemp << endl;
		}
	}
	cout << "outputNum = " << outputNum << endl;
}

int main()
{
	//	PRNG.seed(time(NULL));
	//	string logFile = creatNewFile("log");
	//	PRNG.seed(getFileNum(logFile.c_str()));


	clock_t start = clock();
	clock_t testsec = 0;
	//read_parameter("DNA_system_parameter.txt");
	read_parameter("DNA_system_parameter.txt");
	//K = 30000;

	//n = K + S + H;
	dec_type = "tri";
	set_distri_LT();

	int tot_cases = 10000;
	int	cnt_err_cases = 0;
	cout << "inner code = " << inner_code << endl;
	cout << "run length = " << max_run << endl;




	//stringstream ss;
	//ifstream ifs, ifs2;
	//ofstream ofs;
	//ifs.open("c:\\users\\xinmatrix\\desktop\\data_dec\\1_R10_Symbol.txt", ios::in);
	//ifs2.open("c:\\users\\xinmatrix\\desktop\\data_dec\\1_dec_with_idx.txt", ios::in);
	//unordered_map<int, VI> symbol;
	//int not_exist = 0, exist = 0;
	//long long dis = 0;
	////ofs.open("c:\\users\\xinmatrix\\desktop\\data_dec\\1_dec_with_idx.txt", ios::out);
	//if (!ifs.is_open() || !ifs2.is_open())
	//	assert(0);
	//string stemp1, stemp2;
	//VI vitemp;
	//int itemp, itemp2;
	//int maxline = 11202021;
	//int linenum = 0;
	//int right = 0, wrong = 0;
	//time_t starttime = clock();
	//while (getline(ifs, stemp1)) {
	//	stringstream ss;
	//	ss << stemp1;
	//	ss >> itemp;
	//	vitemp.clear();
	//	for (int i = 0; i < 15; ++i) {
	//		ss >> itemp2;
	//		vitemp.push_back(itemp2);
	//	}
	//	symbol[itemp] = vitemp;
	//}

	//clock_t startTime = clock();
	//while (getline(ifs2, stemp1)) {
	//	++linenum;
	//	if (!(linenum % 1000000)) {
	//		cout << 1.0 * linenum / maxline << ", " << 1.0 * (clock() - startTime) / CLK_TCK;
	//		cout << ", " << 1.0 * dis / exist / l << endl;
	//	}
	//	stringstream ss;
	//	ss << stemp1;
	//	ss >> itemp;
	//	if (!itemp)
	//		continue;
	//	vitemp.clear();
	//	ss >> itemp;
	//	auto iter = symbol.find(itemp);
	//	if (iter == symbol.end()) {
	//		++not_exist;
	//		continue;
	//	}
	//	++exist;
	//	for (int i = 0; i < 15; ++i) {
	//		ss >> itemp2;
	//		vitemp.push_back(itemp2);
	//	}
	//	dis += cal_hanmming_dis(decompress(iter->second, l), decompress(vitemp, l));
	//}
	//cout << "avg dis of undetectable error : " << 1.0 * dis / exist / l << endl;
	//cout << "exist : " << exist << endl;
	//cout << "not exist : " << not_exist << endl;
	//return 0;


	//write_file();
	//return 0;
	//CopyDistributionAfterDecoding2();
	//return 0;
	// testRealCoverage();
	// return 0;
	//SPCC_N_foutain_dec_from_txt();//test the effectiveness of SPCC
	//return 0;
	fountain_dec_from_txt();
	return 0;
	//L_decode_from_txt();
	//return 0;
	//SPCC_dec_info_from_txt();
	//return 0;
	SPCC_N_foutain_dec_from_txt();
	//generate_detect_only_file();
	return 0;







	//bool validSeed = 0;
	//int seed;


	//VD distri_LT;



	//double CON_basis_AP;


	for (int times = 1; times <= tot_cases; ++times) {
		//PRNG.seed(times);
		PRNG.seed(3);
		VVI info = rand_VVI(K, (l + BWI - 1) / BWI);
		inte.clear();
		if (l % BWI != 0) {
			int mask = (1 << l % BWI) - 1;
			for (int i = 0; i < info.size(); ++i) {
				info[i].back() &= mask;//make the tail bits be 0
			}
		}
		VVI write_oligo = write(info);
		continue;
		//		PRNG.seed(time(NULL));
		read_EES = read(write_oligo);

		summary_inner_code(read_EES);

		read_EES = unique_seed(read_EES);

		summary_inner_code(read_EES);

		int max_in_occ = max_incorrect_occurrence(read_EES);

		CN_adj.clear();
		CN_val.clear();

		CN_con.clear();
		CN_err.clear();
		//copy_EES.clear();

		for (int i = 0; i < S + H; ++i) {
			CN_adj.push_back(R10_LH_adj[i]);
			CN_val.push_back(R10_LH_val[i]);
			CN_err.push_back(0);
			CN_con.push_back(1);
			//for (int j = 0; j < read_EES[i].size(); ++j) {
			//	copy_EES[i].push_back(1000);
			//}
		}
		for (int i = 0; i < read_EES.size(); ++i) {
			//			if (read_EES[i].back() > max_in_occ){//correct encoded symbol
			//				CN_adj.push_back(R10_xor_VN(read_EES[i][0]));
			//				CN_val.push_back(VI(read_EES[i].begin() + 1, read_EES[i].end() - 1));
			//			}
			if (read_EES[i].back() > 0) {
				CN_err.push_back(0);
			}
			else if (read_EES[i].back() < 0) {
				CN_err.push_back(1);
			}
			CN_con.push_back(0);
			CN_adj.push_back(R10_xor_VN(read_EES[i][0]));
			CN_val.push_back(VI(read_EES[i].begin() + 1, read_EES[i].end() - 1));
			//for (int j = 0; j < read_EES[i].size(); ++j) {
			//	copy_EES[i + S + H].push_back(read_EES[i][j]);
			//}
		}

		m = CN_adj.size();

		cout << "read_EES.size()=" << read_EES.size() << endl;
		cout << "m=" << m << endl;

		CN_pos_idx.resize(m);
		for (int i = 0; i < m; ++i) {
			CN_pos_idx[i] = i;
		}
		VN_idx_pos.resize(n);
		for (int i = 0; i < VN_idx_pos.size(); ++i) {
			VN_idx_pos[i] = i;
		}

		copy_EES.clear();
		copy_EES.resize(m);//Ҫ�ȸ����������ڴ���ܽ���copy_EES[i]��push_back()����
		for (int i = 0; i < S + H; ++i) {
			for (int j = 0; j < read_EES[0].size() - 1; ++j) {
				copy_EES[i].push_back(1);
			}
			copy_EES[i].push_back(m + 1);
		}
		for (int i = S + H; i < m; ++i) {
			for (int j = 0; j < read_EES[0].size(); ++j) {
				copy_EES[i].push_back(read_EES[i - S - H][j]);
			}
		}

		stable_sort(CN_pos_idx.begin(), CN_pos_idx.end(), cmp1);

		int num = BFA();
		if (num == 1 || num == 2) {
			cnt_err_cases += 1;
		}
		//cnt_err_cases += !decBEC();

		cerr << "cases: " << times << "  FE: " << cnt_err_cases;
		cerr << "  FER: " << 1.0 * cnt_err_cases / times << endl;

		//		FILE *fp = fopen(logFile.c_str(), "w");
		//		fprintf(fp, "inner_code:    %c\n", inner_code);
		//		fprintf(fp, "LT decoding:   %s\n", use_BP ? "BP" : "SGE");
		//		fprintf(fp, "cases:         %d\n", times);
		//		fprintf(fp, "err_cases:     %d\n", cnt_err_cases);
		//		fprintf(fp, "FER:           %.6f\n", 1.0 * cnt_err_cases / times);
		//		fclose(fp);
		//		makeCopy(logFile.c_str());
	}
	testsec = clock() - start;
	cerr << "testsec: " << 1.0 * testsec / CLOCKS_PER_SEC << endl;
	return	0;
}

