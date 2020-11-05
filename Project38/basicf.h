#pragma once
#ifndef BASICF_H_
#define BASICF_H_

#include <vector>
#include <complex>


void flash(const std::string&, size_t, size_t);
std::fstream fileOpener(const std::string&, const std::string&, std::ios_base::openmode);
std::vector<std::string> loadString(const std::string&, const std::string&);
std::vector<size_t> loadInteger(const std::string&, const std::string&);
void setComplexNumericZerotoActualZero(std::complex<double>&);
void roundComplex(std::complex<double>&, int);
std::string convertPhasestoString(const std::vector<std::complex<double>>&, const std::vector<std::complex<double>>& vcd = {});

#endif // !BASICF_H_
