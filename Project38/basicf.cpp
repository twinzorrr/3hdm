#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <Windows.h>
#include "basicf.h"

using boost::lexical_cast;


void flash(const string& text, size_t width, size_t time) {

	cout << setw(width) << left << text << "\r";
	fflush(stdout);
	Sleep(time);

}

void checkFile(fstream& file, const string& filename) {

	if (file.fail()) { cerr << "File " << filename << " Acces Denied!" << endl; exit(EXIT_FAILURE); }

}

fstream fileOpener(const string& path, const string& filename, ios_base::openmode mode) {

	fstream ifile;

	ifile.open(path + filename, mode);
	checkFile(ifile, filename);

	return ifile;

}

vector<string> loadString(const string& path, const string& filename) {

	fstream ifile;
	string line;
	vector<string> vline;

	ifile = fileOpener(path, filename, ios::in);
	while (getline(ifile, line)) vline.push_back(line);
	ifile.close();

	return vline;

}

vector<size_t> loadInteger(const string& path, const string& filename) {

	fstream ifile;
	string line;
	vector<size_t> vn;

	ifile = fileOpener(path, filename, ios::in);
	while (getline(ifile, line)) vn.push_back(lexical_cast<size_t>(line));
	ifile.close();

	return vn;

}

void setComplexNumericZerotoActualZero(complex<double>& cd) {

	if (abs(real(cd)) < 10 * DBL_EPSILON) cd.real(0.0);
	if (abs(imag(cd)) < 10 * DBL_EPSILON) cd.imag(0.0);

}

void roundComplex(complex<double>& cd, int dec) {

	stringstream ss_real, ss_imag;
	ss_real << setprecision(dec) << fixed << real(cd);
	ss_imag << setprecision(dec) << fixed << imag(cd);
	cd.real(stod(ss_real.str())); cd.imag(stod(ss_imag.str()));

}

string convertPhasestoString(vector<complex<double>>& vcd1, vector<complex<double>>& vcd2) {

	string s;
	s = "[ ";
	for (auto& i : vcd1) { setComplexNumericZerotoActualZero(i); roundComplex(i, 6); s += (lexical_cast<string>(i)) + " "; }
	s += "| ";
	for (auto& i : vcd2) { setComplexNumericZerotoActualZero(i); roundComplex(i, 6); s += (lexical_cast<string>(i)) + " "; }
	s += "]";
	return s;

}