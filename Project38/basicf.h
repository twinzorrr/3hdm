#pragma once
#ifndef BASICF_H_
#define BASICF_H_

#include <vector>
#include <complex>

using namespace std;


void flash(const string&, size_t, size_t);
fstream fileOpener(const string&, const string&, ios_base::openmode);
vector<string> loadString(const string&, const string&);
vector<size_t> loadInteger(const string&, const string&);
string convertPhasestoString(vector<complex<double>>&, vector<complex<double>>&);

#endif // !BASICF_H_
