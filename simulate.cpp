#include "simulate.h"


char outer_code;
char inner_code;
int N;
int P;
int deg;
int dmin;
//bool use_BP; defined in decBEC.cpp
double RSD_delta;
double RSD_c;
int n;
int write_size;
int read_size;
int bw_seed;
double c_gc;
int max_run;
bool check_GC_run;
bool detect_only;
double pi, pd, ps;
VI ATCG_to_int(4, 0);
char int_to_ATCG[] = {'A', 'T', 'C', 'G'};
VVD transP(4, VD(4));

int	m;
int	l;
VVI	CN_adj;
VVI VN_adj;
VVI	CN_val;

int SRT_redun() {
	if (max_run == 3)
		return 10;
	else if (max_run == 4)
		return 4;
	else
		assert(0);
}

void set_parameter(const string& par, const string& val) {
	istringstream iss(val);
	if (par == "outer_code") {
		while (iss >> outer_code) {
			if (outer_code == 'R') {
				return;
			}
		}
		assert(false);
	}
	else if (par == "inner_code") {
		while (iss >> inner_code) {
			if (inner_code == 'R' || inner_code == 'L' || inner_code == 'S' || inner_code == 'P') {
				return;
			}
		}
		assert(false);
	}
	else if (par == "N") {
		iss >> N;
	}
	else if (par == "P") {
		iss >> P;
	}
	else if (par == "deg") {
		iss >> deg;
	}
	else if (par == "dmin") {
		iss >> dmin;
	}
	else if (par == "use_BP") {
		iss >> use_BP;
	}
	else if (par == "RSD_delta") {
		iss >> RSD_delta;
	}
	else if (par == "RSD_c") {
		iss >> RSD_c;
	}
	else if (par == "n") {
		if (outer_code == 'R') {
			iss >> K;
			fun();
		}
		else {
			iss >> n;
		}
	}
	else if (par == "write_size") {
		iss >> write_size;
	}
	else if (par == "read_size") {
		iss >> read_size;
	}
	else if (par == "bw_seed") {
		iss >> bw_seed;
	}
	else if (par == "c_gc") {
		iss >> c_gc;
	}
	else if (par == "max_run") {
		iss >> max_run;
	}
	else if (par == "check_GC_run") {
		iss >> check_GC_run;
	}
	else if (par == "detect_only") {
		iss >> detect_only;
	}
	else if (par == "pi") {
		iss >> pi;
	}
	else if (par == "pd") {
		iss >> pd;
	}
	else if (par == "ps") {
		iss >> ps;
	}
	else if (par == "ATCG_to_bits") {
		string tmp = "ATCG";
		for (int i = 0; i < 4; ++i) {
			string bin;
			iss >> bin;
			assert(bin.size() == 2);
			ATCG_to_int[i] = (bin[0] - '0') * 2 + (bin[1] - '0');
			int_to_ATCG[ATCG_to_int[i]] = tmp[i];
		}
		//		output(ATCG_to_int, "ATCG_to_int");
	}
	else if (par == "ATCG_transition_probability") {
		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j) {
				iss >> transP[ATCG_to_int[i]][ATCG_to_int[j]];
			}
		}
		for (int i = 0; i < 4; ++i) {
			normalize(transP[i]);
			for (int j = 1; j < 4; ++j) {
				transP[i][j] += transP[i][j - 1];
			}
		}
		//		output(transP, "transP");
	}
	else {
		cerr << "Unkown parameter: " << par << endl;
	}
}



void read_parameter(const char* fileName) {
	FILE* fp = fopen(fileName, "r");
	if (fp == NULL) {
		cerr << "Cannot open " << fileName << endl;
		assert(false);
	}

	char ch;
	bool new_line = true;
	while ((ch = fgetc(fp)) != EOF) {
		if (ch == '\n') {
			new_line = true;
			continue;
		}

		if (!new_line || ch == '\t' || ch == ' ') {
			continue;
		}

		new_line = false;
		if (ch != '#') {
			continue;
		}

		//read parameter
		//assert(fscanf_s(fp, "%s", data_buff) != EOF);
		//cout << "read" << endl;
		fscanf(fp, "%s", data_buff);
		//cout << "data_buff = " << data_buff << endl;
		int i = 0;
		while (data_buff == NULL) {
			cout << i++;
			fscanf(fp, "s", data_buff);
		}
		//cout << "data_buff = " << data_buff << endl;
		
		//read setting values
		bool value_start = false;
		bool value_end = false;
		string value;
		while ((ch = fgetc(fp)) != EOF) {
			//cout << "ch = " << ch << endl;
			if (ch == '{') {
				value_start = true;
			}
			else if (ch == '}') {
				value_end = true;
				break;
			}
			else if (value_start) {
				value += ch;
			}
			else if (ch != '\t' && ch != ' ' && ch != '\n') {
				cerr << "Non-blank charachter occurs between par and {" << endl;
				cout << "ch = " << ch << endl;
				assert(false);
			}
		}

		if (!value_start || !value_end) {
			cerr << "Have not found { or }" << endl;
			assert(false);
		}
		else {
			set_parameter(data_buff, value);
		}
	}

	assert(bw_seed <= 32 && (1LL << bw_seed) >= write_size);

	l = R10_info_length() - SRT_redun();
	//l = R10_info_length();
	cout << "bit-width of R10 source symbol l: " << l << endl << endl;

	fclose(fp);
}



int max_run_binary(const VB& A) {
	assert(!A.empty());
	int run = 1, res = 1;
	for (int i = 1; i < A.size(); ++i) {
		if (A[i] == A[i - 1]) {
			++run;
			if (run > res) {
				res = run;
			}
		}
		else {
			run = 1;
		}
	}
	return	res;
}

bool valid_GC_run(const VB& oligo) {
	assert(oligo.size() % 2 == 0);
	int run = 1;
	int nt = oligo[0] * 2 + oligo[1], pre_nt = nt;
	int cnt_GC = (ATCG_to_int[2] == nt || ATCG_to_int[3] == nt);
	for (int i = 2; i < oligo.size(); i += 2) {
		nt = oligo[i] * 2 + oligo[i + 1];
		cnt_GC += (ATCG_to_int[2] == nt || ATCG_to_int[3] == nt);

		if (nt == pre_nt) {
			++run;
			if (run > max_run) {
				//cout << "i == " << i << ", " << run << ", " << nt << endl;
				return	false;
			}
		}
		else {
			run = 1;
			pre_nt = nt;
		}
	}

	return	abs(2.0 * cnt_GC / oligo.size() - 0.5) <= c_gc;
}

bool valid_run_int(int seed)
{
	int run = 1;
	int nt = seed & 3, pre_nt = nt;
	for (int i = 2; i < bw_seed; i += 2) {
		seed >>= 2;
		nt = seed & 3;
		if (nt == pre_nt) {
			++run;
			if (run > max_run) {
				return	false;
			}
		}
		else {
			run = 1;
			pre_nt = nt;
		}
	}

	return	true;
}

VB interleave(const VB& A, const VB& B) {
	assert(A.size() == B.size());
	VB res(A.size() * 2);
	for (int i = 0; i < A.size(); ++i) {
		res[i * 2] = A[i];
		res[i * 2 + 1] = B[i];
	}
	return	res;
}

VVB de_interleave(const VB& A) {
	assert(A.size() % 2 == 0);
	VVB res(2, VB(A.size() / 2));
	for (int i = 0; i < A.size(); ++i) {
		res[i % 2][i / 2] = A[i];
	}
	return	res;
}

VI	compress_ES(const VB& A) {
	VI res = compress(VB(A.begin(), A.begin() + bw_seed));
	VI tmp = compress(VB(A.begin() + bw_seed, A.end()));
	res.insert(res.end(), tmp.begin(), tmp.end());
	return	res;
}


VB	decompress_ES(const VI& A) {
	VB res = decompress(VI(A.begin(), A.begin() + 1), bw_seed);
	VB tmp = decompress(VI(A.begin() + 1, A.end()), l);
	res.insert(res.end(), tmp.begin(), tmp.end());
	return	res;
}

VI bin2qua(const VB& A) {
	VVB temp = de_interleave(A);
	VI qua;
	for (int i = 0; i < temp[0].size(); ++i) {
		if (temp[0][i] == 1 && temp[1][i] == 1)
			qua.push_back(3);
		else if (temp[0][i] == 1 && temp[1][i] == 0)
			qua.push_back(2);
		else if (temp[0][i] == 0 && temp[1][i] == 1)
			qua.push_back(1);
		else
			qua.push_back(0);
	}
	return qua;
}

VB qua2bin(const VI& A) {
	VVB temp(2, VB(0));
	for (int i = 0; i < A.size(); ++i) {
		if (A[i] == 3) {
			temp[0].push_back(1);
			temp[1].push_back(1);
		}
		else if (A[i] == 2) {
			temp[0].push_back(1);
			temp[1].push_back(0);
		}
		else if (A[i] == 1) {
			temp[0].push_back(0);
			temp[1].push_back(1);
		}
		else {
			temp[0].push_back(0);
			temp[1].push_back(0);
		}
	}
	return interleave(temp[0], temp[1]);
}

VB string2bin(const string& s)
{
	VB rst;
	for (int i = 0; i < s.size(); ++i) {
		if (s[i] == 'A') {
			rst.push_back(0);
			rst.push_back(0);
		}
		else if (s[i] == 'T') {
			rst.push_back(0);
			rst.push_back(1);
		}
		else if (s[i] == 'C') {
			rst.push_back(1);
			rst.push_back(0);
		}
		else {
			rst.push_back(1);
			rst.push_back(1);
		}
	}

	return rst;
}

string bin2string(const VB& v)
{
	string rst;
	VVB tmp = de_interleave(v);
	for (int i = 0; i < tmp[0].size(); ++i) {
		if (tmp[0][i] == 1 && tmp[1][i] == 1)
			rst += "G";
		else if (tmp[0][i] == 1 && tmp[1][i] == 0)
			rst += "C";
		else if (tmp[0][i] == 0 && tmp[1][i] == 1)
			rst += "T";
		else
			rst += "A";
	}

	return rst;
}

VVI write(const VVI& info) {
	ofstream ofs;
	//ofs.open("C:\\Users\\xinmatrix\\Desktop\\Data_dec\\3_R10_Symbol.txt", ios::out);
	//if (!ofs.is_open())
		//assert(0);


	assert(!info.empty());
	int codingScheme = 3;//0, 1, 2, 3 represent L+3, L+4, P+3, P+4 respectively
	VVI res;
	int countEnc = 0;
	while (res.size() < write_size) {
		++countEnc;


		//VB code = decompress_ES(R10_encode(info));
		VI R10_symbol = R10_encode(info);
		VB code = decompress_ES(R10_symbol);


		VB codeInfo(code.begin() + bw_seed, code.end());

		//SRT encoding
		codeInfo = SRT_encode(codeInfo);
		//code = SRT_encode(code);
		//assert(valid_GC_run(codeInfo));
		//if(!valid_GC_run(codeInfo)) {
		//	output(bin2qua(codeInfo), "code quantury");
		//}

		code.resize(bw_seed);
		reverse(code.begin(), code.end());
		code.insert(code.end(), codeInfo.begin(), codeInfo.end());

		string cs = bin2string(code);
		int gc = 0;
		for (auto a : cs) {
			if (a == 'G' || a == 'C')
				++gc;
		}
		int a = gc;


		//code.insert(code.begin(), codeInfo.begin(), codeInfo.end());

		//ECC encoding
		if (inner_code == 'R') {
			code = RS_encode(code);
		}
		else if (inner_code == 'L') {
			
			
			code = L_encode(code);
		}
		else if (inner_code == 'S') {
			code = SVT_encode(code);
			
		}
		else if (inner_code == 'P') {
			code = SPCC_encode(code);

			string cs = bin2string(code);
			int gc = 0;
			for (auto a : cs) {
				if (a == 'G' || a == 'C')
					++gc;
			}
			int a = gc;
		}
		if (!code.empty()) {
			code.push_back((codingScheme >> 1)& 1);
			code.push_back(codingScheme & 1);
			code.emplace(code.begin(), codingScheme & 1);
			code.emplace(code.begin(), (codingScheme >> 1) & 1);
			
		}

		if (!code.empty() && (!check_GC_run || valid_GC_run(code))) {
			//output(bin2qua(code), "code");
			res.push_back(compress(code));
			//for (int i = 0; i < R10_symbol.size(); ++i)
			//	ofs << R10_symbol[i] << " ";
			//ofs << endl;
		}
	}
	cout << "count encoding = " << countEnc << endl;
	

	//FILE* fp = fopen("output\\_C4.txt", "w");
	//if (fp != NULL) {
	//	for (int i = 0; i < res.size(); ++i) {
	//		VI temp = bin2qua(decompress(res[i], 2 * (N + 2) ));
	//		for (auto a : temp) {
	//			fprintf(fp, "%c", int_to_ATCG[a]);
	//		}
	//		fprintf(fp, "\n");
	//	}
	//fclose(fp);
	//}
	

	//fp = fopen("output\\L+3+0.05info.txt", "w");
	//if (fp != NULL) {
	//	for (int i = 0; i < info.size(); ++i) {
	//		VB temp = decompress(info[i], l);
	//		
	//		for (int j = 0; j < temp.size(); ++j) {
	//			if (temp[j]) fprintf(fp, "1 ");
	//			else fprintf(fp, "0 ");
	//		}
	//		fprintf( fp, "\n\n");
	//	}
	//}
	//fclose(fp);

	return	res;
}

VB add_noise(const VB& oligo) {//oligo in binary form: two binary codewords interleaved
	assert(oligo.size() == 2 * N);

	//	output(de_interleave(oligo), "before noise");

	VB res;
	for (int i = 0; i < oligo.size(); i += 2) {
		int cur_nt = oligo[i] * 2 + oligo[i + 1];
		double p = randu();
		if (p < pi) {//insert an nt before cur_nt
			res.push_back(randb());
			res.push_back(randb());
			res.push_back(oligo[i]);
			res.push_back(oligo[i + 1]);
		}
		else if (p < pi + pd) {//delete cur_nt
			;
		}
		else if (p < pi + pd + ps) {//substitute cur_nt
			p = randu();
			for (int new_nt = 0; new_nt < 4; ++new_nt) {
				if (p < transP[cur_nt][new_nt]) {
					res.push_back(new_nt >> 1);
					res.push_back(new_nt & 1);
					break;
				}
			}
		}
		else {//remain cur_nt
			res.push_back(oligo[i]);
			res.push_back(oligo[i + 1]);
		}
	}

	double p = randu();
	if (p < pi) {//insert an nt at the tail
		res.push_back(randb());
		res.push_back(randb());
	}

	//	output(de_interleave(res), "after noise");
	return	res;
}

VVI read(const VVI& write_oligo) {
	VVI read_EES;
	if (inner_code == 'R') {
		read_EES = RS_read(write_oligo);
	}
	else if (inner_code == 'L') {
		read_EES = L_read(write_oligo);
	}
	else if (inner_code == 'S') {
		read_EES = SVT_read(write_oligo);
	}
	return	read_EES;
}

VVI unique_seed(const VVI& read_EES) {
	VPII pr(read_EES.size());
	for (int i = 0; i < pr.size(); ++i) {
		pr[i].first = read_EES[i][0];//seed
		pr[i].second = i;
	}
	sort(pr.begin(), pr.end());//sort based on seeds

	VVI res;
	for (int i = 0; i < pr.size();) {
		int occur_max = abs(read_EES[pr[i].second].back());
		int idx = i;
		int cnt = 1;

		int j = i;
		while (++j < pr.size() && read_EES[pr[j].second][0] == read_EES[pr[i].second][0]) {
			if (abs(read_EES[pr[j].second].back()) > occur_max) {
				occur_max = abs(read_EES[pr[j].second].back());
				idx = j;
				cnt = 1;
			}
			else if (abs(read_EES[pr[j].second].back()) == occur_max) {//oligo���������ﱻ������
				++cnt;
			}
		}

		if (cnt == 1) {//the maximum copy is unique
			res.push_back(read_EES[pr[idx].second]);
		}
		i = j;
	}

	return	res;
}

int	max_incorrect_occurrence(const VVI& read_EES) {
	int res = 0;
	for (int i = 0; i < read_EES.size(); ++i) {
		res = min(res, read_EES[i].back());
	}
	return	-res;
}

void summary_inner_code(const VVI& read_EES) {
	int cnt_valid = 0;//number of valid output oligos
	int exist = 0;//number of orignal encoded symbols that appear at output
	int cnt_incorrect = 0;//number of incorrect copies
	int cnt_correct = 0;//number of correct copies
	int max_in_occ = max_incorrect_occurrence(read_EES);
	VI	cnt_incorrect_each(max_in_occ + 1, 0);
	VI  cnt_correct_each(max_in_occ + 5, 0);

	for (int i = 0; i < read_EES.size(); ++i) {
		assert(read_EES[i].back() != 0);

		cnt_valid += abs(read_EES[i].back());
		if (read_EES[i].back() < 0) {//incorrect output
			cnt_incorrect -= read_EES[i].back();
			++cnt_incorrect_each[-read_EES[i].back()];
		}
		else {
			cnt_correct += read_EES[i].back();
			++exist;
			if (read_EES[i].back() < cnt_correct_each.size()) {
				++cnt_correct_each[read_EES[i].back()];
			}
		}
	}
	cnt_correct_each[0] = write_size - exist;

	if (inner_code == 'R') {
		cerr << "RS code" << endl << endl;
	}
	else if (inner_code == 'L') {
		cerr << "Levenshtein code" << endl << endl;
	}
	else if (inner_code == 'S') {
		cerr << "SVT code" << endl << endl;
	}
	else {
		cerr << "inner code is incorrect." << endl;
		assert(false);
	}

	cout << "n                 \t" << n << endl;
	cout << "write_size        \t" << write_size << endl;
	cout << "read_size         \t" << read_size << endl;
	cout << "cnt_valid         \t" << cnt_valid << " \t" << 1.0 * cnt_valid / read_size << endl;
	cout << "cnt_incorrect     \t" << cnt_incorrect << " \t" << 1.0 * cnt_incorrect / read_size << endl;
	cout << "cnt_correct       \t" << cnt_correct << " \t" << 1.0 * cnt_correct / read_size << endl;

	cout << endl << "incorrect" << endl;
	cout << "occurrence \tnumber \trate (number/n)" << endl;
	for (int j = 0; j < cnt_incorrect_each.size(); ++j) {
		printf("%d\t\t%d\t%.4f\n", j, cnt_incorrect_each[j], 1.0 * cnt_incorrect_each[j] / n);
	}
	cout << "----" << endl << endl;

	cout << endl << "correct" << endl;
	cout << "occurrence \tnumber \trate \trate_sum \t(write_size - sum_num)/n" << endl;
	double sum = 0.0;
	for (int j = 0; j < cnt_correct_each.size(); ++j) {
		sum += cnt_correct_each[j];
		printf("%d\t\t%d\t%.4f\t%.4f\t\t%.4f\n", j, cnt_correct_each[j], \
			1.0 * cnt_correct_each[j] / n, sum / n, (write_size - sum) / n);
		if (j + 1 == cnt_incorrect_each.size()) {
			cout << endl;
		}
	}
	cout << endl << "--------------------------------------" << endl << endl;
}

void SRT_enc_4A_block(VI& diff_info) {
	assert(max_run >= 2 && diff_info.size());
	assert(diff_info.size() <= pow(4, max_run - 1) * 3 + max_run - 1);
	diff_info.push_back(0);
	if (diff_info.size() <= max_run)
		return;
	
	int pos = 1;
	int sz = 1;
	while (pos < diff_info.size()) {
		
		int idx = pos - 1;
		while (++idx < diff_info.size() && diff_info[idx] == 0 && idx < pos + max_run);
		if (idx == pos + max_run) {
			int addr = sz - 1;
			for (int i = 0; i < max_run; ++i) {
				diff_info.push_back(addr & 3);
				addr >>= 2;
			}
			diff_info.back() += 1;
			pos = idx;
		}
		else {
			do {
				diff_info[sz++] = diff_info[pos];
			} while (++pos < idx);
		}
		
		//take example of rll = 3
		//pos <-> c1c2c3 (c1 = least significant digit)
		//0   <-> 1 0 0
		//1   <-> 1 0 1
		//2   <-> 1 0 2
		// ...
		//47  <-> 3 3 3
		//so, pos + 16 -> c1c2c3 
		
		//for (auto a : diff_info)
		//	cout << a << "  ";
		//cout << endl;
	}
	//return prefix;
	diff_info.resize(sz);
}

void SRT_dec_4A_block(VI& diff_info)
{
	if (diff_info.empty())
		return;
	while (diff_info.back() != 0) {
		int addr = 0;
		diff_info.back() -= 1;
		for (int i = 0; i < max_run; ++i) {
			addr += diff_info.back() << 2 * (max_run - i - 1);
			//addr += diff_info.back() << 2 * i;
			diff_info.pop_back();
		}
		//assert(addr < diff_info.size());
		if (addr >= diff_info.size()) {
			diff_info.resize(0);
			return;
		}
		for (int i = 0; i < max_run; ++i) {
			diff_info.emplace(diff_info.begin() + addr + 1, 0);
		}
	}
	diff_info.pop_back();
}



VB SRT_encode(VB const& infoBin) {
	VI info = bin2qua(infoBin);
	
	int blockLen = pow(4, max_run - 1) * 3 + max_run - 1;
	VVI info2;
	VI temp;
	int t = 0;
	while (t < info.size() && t < blockLen) {
		temp.push_back(info[t++]);
	}
	for (int i = temp.size() - 1; i > 0; --i)
		temp[i] ^= temp[i - 1];
	SRT_enc_4A_block(temp);
	for (int j = 1; j < temp.size(); ++j)
		temp[j] ^= temp[j - 1];
	info2.push_back(temp);
	while (t < info.size()) {
		temp.resize(0);
		while ((int)temp.size() < blockLen - max_run && t < info.size())
			temp.push_back(info[t++]);
		
		//temp.insert(temp.begin(), info2.back().end() - max_run, info2.back().end());
		//for (int i = 0; i < max_run; ++i) {
		//	info2.back().pop_back();
		//}
		for (int i = 0; i < max_run; ++i) {
			temp.emplace(temp.begin(), info2.back().back());
			info2.back().pop_back();
		}

		for (int i = temp.size() - 1; i > 0; --i)
			temp[i] ^= temp[i - 1];
		SRT_enc_4A_block(temp);
		for (int j = 1; j < temp.size(); ++j)
			temp[j] ^= temp[j - 1];
		info2.push_back(temp);
	}
	//for (auto a : info2)
	//	for (auto b : a)
	//		cout << b << ", ";
	//cout << endl;
	temp.resize(0);
	for (auto a : info2)
		for (auto b : a)
			temp.push_back(b);
	//cout << "info size = " << info.size() << endl;
	//cout << "temp size = " << temp.size() << endl;
	return qua2bin(temp);

	//for (int i = 0; i < blockNum; ++i) {
	//	temp.resize(0);
	//	for (int j = i * blockLen; j < (i + 1) * blockLen && j < info.size(); ++j)
	//		temp.push_back(info[j]);
	//	info2.push_back((temp));
	//}
	//for (int i = 0; i < info2.size(); ++i)
	//	for (int j = 0; j < info2[i].size() - 1; ++j) 
	//		info2[i][j] = info2[i][j + 1] - info2[i][j];
	//for (auto a : info2)
	//	a = SRT_enc_4A_block(a);
	//for (int i = 0; i < info2.size(); ++i)
	//	for (int j = info2[i].size() - 2; j >= 0; --j)
	//		info2[i][j] = info2[i][j + 1] + info2[i][j];
	//temp.resize(0);
	//for (auto a : info2)
	//	for (auto b : a)
	//		temp.push_back(b);
	//return temp;
}

VB SRT_decode(VB const& cwBin)
{
	VI cw = bin2qua(cwBin);
	
	int blockLen = pow(4, max_run - 1) * 3 + max_run;
	int sz = 0;
	if (cw.size() <= blockLen)
		sz = cw.size();
	else {
		sz = cw.size() - blockLen;
		while (sz > blockLen - max_run)
			sz -= blockLen - max_run;
		sz += max_run;
	}

	int t = cw.size() - 1;
	VI temp;
	VVI res;
	while (t >= 0 && temp.size() < sz) {
		temp.emplace(temp.begin(), cw[t--]);
	}
	for (int i = temp.size() - 1; i > 0; --i)
		temp[i] ^= temp[i - 1];
	SRT_dec_4A_block(temp);
	if (temp.empty())
		return VB(0);
	for (int j = 1; j < temp.size(); ++j)
		temp[j] ^= temp[j - 1];
	res.emplace(res.begin(), temp);
	while (t >= 0) {
		temp.resize(0);
		while (t >= 0 && (int) temp.size() < blockLen - max_run) {
			temp.emplace(temp.begin(), cw[t--]);
		}
		for (int i = 0; i < max_run; ++i) {
			temp.push_back(res.front().front());
			res.front().erase(res.front().begin());
		}
		for (int i = temp.size() - 1; i > 0; --i)
			temp[i] ^= temp[i - 1];
		SRT_dec_4A_block(temp);
		if (temp.empty())
			return VB(0);
		for (int j = 1; j < temp.size(); ++j)
			temp[j] ^= temp[j - 1];
		res.emplace(res.begin(), temp);
	}
	temp.resize(0);
	for (auto a : res)
		for (auto b : a)
			temp.push_back(b);
	return qua2bin(temp);
}

