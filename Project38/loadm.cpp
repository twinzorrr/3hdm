#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <regex>
#include "loadm.h"
#include "loadss.h"

using boost::lexical_cast;

static complex<double> i(0.0, 1.0);
static const double pi = acos(-1.0);


regex re() {

	string s1, s2, s3, s4, s5, s6, s7;
	s1 = "\\s?-?(\\d+)/(\\d+)\\*E\\((\\d+)\\)\\^(\\d+)";
	s2 = "\\s?-?(\\d+)/(\\d+)\\*E\\((\\d+)\\)";
	s3 = "\\s?-?(\\d+)\\*E\\((\\d+)\\)\\^(\\d+)";
	s4 = "\\s?-?(\\d+)\\*E\\((\\d+)\\)";
	s5 = "\\s-?(\\d+)";
	s6 = ".?E\\((\\d+)\\)\\^(\\d+)";
	s7 = ".?E\\((\\d+)\\)";
	regex s(s1 + "|" + s2 + "|" + s3 + "|" + s4 + "|" + s5 + "|" + s6 + "|" + s7);

	return s;

}

vector<double> d(const smatch& res) {

	int i;
	double d[4];

	i = 0; d[1] = 1.0; d[2] = 1.0; d[3] = 1.0;
	for (size_t r = 1; r < res.size(); r++) {
		if (res[r] != "") { d[i] = lexical_cast<double>(res[r]); i++; }
	}

	return { d[0],d[1],d[2],d[3] };
}

complex<double> e(const string& s, double d1, double d2, double d3, double d4) {

	complex<double> c, d;
	bool con;

	d = { 2.0 * d2 / d1, 0.0 }; c = { 1.0, 0.0 };
	if (s.find("/") != string::npos) { d = { 2.0 * d4 / d3, 0.0 }; c = { d1 / d2, 0.0 }; }
	else if (s.find("*") != string::npos) { d = { 2.0 * d3 / d2, 0.0 }; c = { d1, 0.0 }; }

	con = (s[0] == '-') | (s[1] == '-');
	if (con) return -1.0 * c * exp(i * pi * d);
	return c * exp(i * pi * d);

}

vector<MatrixXcd> loadM(fstream& file, int dim, int prec) {

	regex r(re());
	string::const_iterator start;
	string line, sec;
	int k;
	complex<double> f, g(0.0, 0.0);
	stringstream ss;
	MatrixXcd a;
	vector<MatrixXcd> va;

	k = 0;

	while (getline(file, line)) {
		smatch res;
		start = line.cbegin();
		while (regex_search(start, line.cend(), res, r)) {
			k++;
			start += res.position() + res.length();
			sec = res.str(0);
			if (sec.find("E") != string::npos) {
				ostringstream oss;
				f = e(sec, d(res)[0], d(res)[1], d(res)[2], d(res)[3]);
				if ((*start == '-') | (*start == '+')) { k--; g += f; continue; }
				f += g;
				oss << setprecision(prec) << f;
				sec = " " + oss.str();
				g = 0.0;
			}
			if (k == dim * dim) { ss << sec; loadMatrix(ss, a, dim); va.push_back(a); k = 0; ss.str(""); ss.clear(); continue; }
			sec = sec;
			ss << sec;
		}
	}

	return va;

}