#pragma once
#ifndef LOADM_H_
#define LOADM_H

#include "loadss.h"

template<typename T> using Real = typename Eigen::NumTraits<T>::Real;


namespace load {

	template<typename _scalar>
	std::vector<Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>> loadM(std::fstream& file, int dim, int prec) {

		assert(Eigen::NumTraits<_scalar>::IsComplex);

		std::regex r(detail::re());
		std::string::const_iterator start;
		std::string line, sec;
		int k = 0;
		std::vector<Real<_scalar>> vsc;
		_scalar f, g((Real<_scalar>)0, (Real<_scalar>)0);
		std::stringstream ss;
		Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic> a;
		std::vector<Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>> va;

		while (getline(file, line)) {
			std::smatch res;
			start = line.cbegin();
			while (std::regex_search(start, line.cend(), res, r)) {
				k++;
				start += res.position() + res.length();
				sec = res.str(0);
				if (sec.find("E") != std::string::npos) {
					std::ostringstream oss;
					vsc = detail::d<Real<_scalar>>(res);
					f = detail::e<_scalar>(sec, vsc[0], vsc[1], vsc[2], vsc[3]);
					if ((*start == '-') || (*start == '+')) { k--; g += f; continue; }
					f += g;
					oss << std::setprecision(prec) << f;
					sec = " " + oss.str();
					g = { (Real<_scalar>)0, (Real<_scalar>)0 };
				}
				if (k == dim * dim) { ss << sec; detail::loadMatrix(ss, a, dim); va.push_back(a); k = 0; ss.str(""); ss.clear(); continue; }
				sec = sec;
				ss << sec;
			}
		}

		return va;

	}

}

#endif // !LOADM_H_