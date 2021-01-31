#pragma once
#ifndef BASICF_H_
#define BASICF_H_

#include "Eigen/Dense"

template<typename T> using Real = typename Eigen::NumTraits<T>::Real;


namespace basic {

	void flash(const std::string&, size_t, size_t);
	std::fstream fileOpener(const std::string&, const std::string&, std::ios_base::openmode);
	std::vector<std::string> loadString(const std::string&, const std::string&);
	std::vector<size_t> loadInteger(const std::string&, const std::string&);
	template<typename _scalar> void setComplexNumericZerotoActualZero(_scalar&);
	template<typename _scalar> bool isNotUniqueComplex(const _scalar, const _scalar);
	template<typename _scalar> void updateGlobalPhaseStandardization(std::vector<_scalar>&, _scalar);
	template<typename _scalar> void setGlobalPhaseStandardization(std::vector<_scalar>&, std::vector<_scalar>&);
	template<typename _scalar> const std::string convertPhasestoString(std::vector<_scalar>&, std::vector<_scalar>&);
	template<typename _scalar> const std::string convertPhasestoString(std::vector<_scalar>&);

	template<typename _scalar>
	void setComplexNumericZerotoActualZero(_scalar& sc) {

		if (abs(real(sc)) < (Real<_scalar>)0.000001) sc.real((Real<_scalar>)0);
		if (abs(imag(sc)) < (Real<_scalar>)0.000001) sc.imag((Real<_scalar>)0);

	}

	template<typename _scalar>
	bool isNotUniqueComplex(const _scalar sc1, const _scalar sc2) {

		return abs(sc1.real() - sc2.real()) < (Real<_scalar>)0.00001 && abs(sc1.imag() - sc2.imag()) < (Real<_scalar>)0.00001;

	}

	template<typename _scalar>
	void updateGlobalPhaseStandardization(std::vector<_scalar>& vscst, _scalar sc) {

		setComplexNumericZerotoActualZero(sc);
		for (const auto& i : vscst) {
			if (isNotUniqueComplex(i, sc)) return;
		}
		vscst.push_back(sc);

	}

	template<typename _scalar>
	void setGlobalPhaseStandardization(std::vector<_scalar>& vscst, std::vector<_scalar>& vsc) {

		for (const auto& i : vsc) updateGlobalPhaseStandardization(vscst, i);
		for (const auto& i : vscst) std::replace_if(vsc.begin(), vsc.end(), [=](_scalar sc) { return basic::isNotUniqueComplex(i, sc); }, i);

	}

	template<typename _scalar>
	const std::string convertPhasestoString(std::vector<_scalar>& vsc1, std::vector<_scalar>& vsc2) {

		using boost::lexical_cast;

		static std::vector<_scalar> vscst;
		basic::setGlobalPhaseStandardization(vscst, vsc1);

		std::string s;
		s = "[ ";
		for (const auto& i : vsc1) { s += (lexical_cast<std::string>(i)) + " "; }
		if (vsc2.size()) {
			basic::setGlobalPhaseStandardization(vscst, vsc2);
			s += "| ";
			for (const auto& i : vsc2) { s += (lexical_cast<std::string>(i)) + " "; }
		}
		s += "]";
		return s;

	}

	template<typename _scalar>
	const std::string convertPhasestoString(std::vector<_scalar>& vsc1) {

		std::vector<_scalar> vsc2 = {};
		return convertPhasestoString(vsc1, vsc2);

	}

}

#endif // !BASICF_H_