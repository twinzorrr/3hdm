#include <iostream>
#include <iomanip>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <Windows.h>
#include "basicf.h"


namespace basic {

	void flash(const std::string& text, size_t width, size_t time) {

		std::cout << std::setw(width) << std::left << text << "\r";
		fflush(stdout);
		Sleep(time);

	}

	void checkFile(std::fstream& file, const std::string& filename) {

		if (file.fail()) { std::cerr << "File " << filename << " Access Denied!" << std::endl; exit(EXIT_FAILURE); }

	}

	std::fstream fileOpener(const std::string& path, const std::string& filename, std::ios_base::openmode mode) {

		std::fstream ifile;

		ifile.open(path + filename, mode);
		checkFile(ifile, filename);

		return ifile;

	}

	std::vector<std::string> loadString(const std::string& path, const std::string& filename) {

		std::fstream ifile;
		std::string line;
		std::vector<std::string> vline;

		ifile = fileOpener(path, filename, std::ios::in);
		while (getline(ifile, line)) vline.push_back(line);
		ifile.close();

		return vline;

	}

	std::vector<size_t> loadInteger(const std::string& path, const std::string& filename) {

		using boost::lexical_cast;

		std::fstream ifile;
		std::string line;
		std::vector<size_t> vn;

		ifile = fileOpener(path, filename, std::ios::in);
		while (getline(ifile, line)) vn.push_back(lexical_cast<size_t>(line));
		ifile.close();

		return vn;

	}

	void setComplexNumericZerotoActualZero(std::complex<double>& cd) {

		if (abs(real(cd)) < 10 * DBL_EPSILON) cd.real(0.0);
		if (abs(imag(cd)) < 10 * DBL_EPSILON) cd.imag(0.0);

	}

	void roundComplex(std::complex<double>& cd, int dec) {

		std::stringstream ss_real, ss_imag;
		ss_real << std::setprecision(dec) << std::fixed << real(cd);
		ss_imag << std::setprecision(dec) << std::fixed << imag(cd);
		cd.real(stod(ss_real.str())); cd.imag(stod(ss_imag.str()));

	}

	std::string convertPhasestoString(const std::vector<std::complex<double>>& vcd1, const std::vector<std::complex<double>>& vcd2) {

		using boost::lexical_cast;

		std::string s;
		s = "[ ";
		for (auto i : vcd1) { setComplexNumericZerotoActualZero(i); roundComplex(i, 6); s += (lexical_cast<std::string>(i)) + " "; }
		if (vcd2.size()) {
			s += "| ";
			for (auto i : vcd2) { setComplexNumericZerotoActualZero(i); roundComplex(i, 6); s += (lexical_cast<std::string>(i)) + " "; }
		}
		s += "]";
		return s;

	}

}