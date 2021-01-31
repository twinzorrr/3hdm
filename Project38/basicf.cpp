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
		std::vector<std::string> vs;

		ifile = fileOpener(path, filename, std::ios::in);
		while (getline(ifile, line)) vs.push_back(line);
		ifile.close();

		return vs;

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

}