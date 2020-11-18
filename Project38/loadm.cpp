#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/lexical_cast.hpp>
#include <regex>
#include "loadm.h"
#include "loadss.h"

using boost::lexical_cast;

static std::complex<double> i(0.0, 1.0);
static const double pi = acos(-1.0);


std::regex re() {

	std::string s1, s2, s3, s4, s5, s6, s7;
	s1 = "\\s?-?(\\d+)/(\\d+)\\*E\\((\\d+)\\)\\^(\\d+)";
	s2 = "\\s?-?(\\d+)/(\\d+)\\*E\\((\\d+)\\)";
	s3 = "\\s?-?(\\d+)\\*E\\((\\d+)\\)\\^(\\d+)";
	s4 = "\\s?-?(\\d+)\\*E\\((\\d+)\\)";
	s5 = "\\s-?(\\d+)";
	s6 = ".?E\\((\\d+)\\)\\^(\\d+)";
	s7 = ".?E\\((\\d+)\\)";
	std::regex s(s1 + "|" + s2 + "|" + s3 + "|" + s4 + "|" + s5 + "|" + s6 + "|" + s7);

	return s;

}

std::vector<double> d(const std::smatch& res) {

	int i;
	double d[4];

	i = 0; d[1] = 1.0; d[2] = 1.0; d[3] = 1.0;
	for (size_t r = 1; r < res.size(); r++) {
		if (res[r] != "") { d[i] = lexical_cast<double>(res[r]); i++; }
	}

	return { d[0],d[1],d[2],d[3] };
}

std::complex<double> e(const std::string& s, double d1, double d2, double d3, double d4) {

	std::complex<double> c, d;
	bool con;

	d = { 2.0 * d2 / d1, 0.0 }; c = { 1.0, 0.0 };
	if (s.find("/") != std::string::npos) { d = { 2.0 * d4 / d3, 0.0 }; c = { d1 / d2, 0.0 }; }
	else if (s.find("*") != std::string::npos) { d = { 2.0 * d3 / d2, 0.0 }; c = { d1, 0.0 }; }

	con = (s[0] == '-') || (s[1] == '-');
	if (con) return -1.0 * c * exp(i * pi * d);
	return c * exp(i * pi * d);

}

std::vector<Eigen::MatrixXcd> loadM(std::fstream& file, int dim, int prec) {

	std::regex r(re());
	std::string::const_iterator start;
	std::string line, sec;
	int k;
	std::complex<double> f, g(0.0, 0.0);
	std::stringstream ss;
	Eigen::MatrixXcd a;
	std::vector<Eigen::MatrixXcd> va;

	k = 0;

	while (getline(file, line)) {
		std::smatch res;
		start = line.cbegin();
		while (std::regex_search(start, line.cend(), res, r)) {
			k++;
			start += res.position() + res.length();
			sec = res.str(0);
			if (sec.find("E") != std::string::npos) {
				std::ostringstream oss;
				f = e(sec, d(res)[0], d(res)[1], d(res)[2], d(res)[3]);
				if ((*start == '-') || (*start == '+')) { k--; g += f; continue; }
				f += g;
				oss << std::setprecision(prec) << f;
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