#include <iostream>
#include <fstream>
#include "basicf.h"
#include "myvector.h"
#include "mymatrix.h"
#include "group.h"
#include "yukawa.h"


std::vector<Eigen::Matrix3cd> Yukawa::convertSolutionstoMassMatrices() const {

	Eigen::Matrix3cd m;
	std::vector<Eigen::Matrix3cd> vm;

	for (size_t i = 0; i < vv_.size(); i++) {
		m = vv_[i].getMassMatrix() * vv_[i].getMassMatrix().adjoint();
		vm.push_back(m);
	}

	return vm;

}

void Yukawa::printMassRatio(std::vector<double>& vd, std::fstream& file) const {

	file << "Muon electron mass ratio" << std::endl << " min: " << vd[0] << " max: " << vd[1] << std::endl;
	file << "Tau electron mass ratio" << std::endl << " min: " << vd[2] << " max: " << vd[3] << std::endl;
	file << "Tau muon mass ratio" << std::endl << " min: " << vd[4] << " max: " << vd[5] << std::endl;
	file << std::endl << std::endl;

}

Yukawa::Yukawa() : solution_(yo::UNDEFINED) {}

Yukawa::Yukawa(yo solution, std::vector<MyVector>& vmv, const std::vector<std::string>& vs)
	: solution_(solution), vs_(vs), vv_(vmv)
{
	for_each(vmv.begin(), vmv.end(), [](MyVector& mv) { assert(mv.getNumberofElements() == 27); });
}

void Yukawa::setNumericZerotoActualZero() { for (auto& i : vv_) i.setNumericZerotoActualZero(); }

void Yukawa::setFirstNonZeroElementto1() { for (auto& i : vv_) i.setFirstNonZeroElementto1(); }

bool Yukawa::errorGuard(const Group& g) const {

	for (auto i : vv_) {
		if (!g.isEigenvector1(i)) { std::cerr << "Yukawa solution is not a eigenvector!" << std::endl; return true; }
		if (i.isApprox(MyVector(Eigen::VectorXcd::Zero(27)), 0.01)) { std::cerr << "Yukawa solution is a zero vector!"; return true; }
	}

	for (size_t i = 0; i < getSize() - 1; i++) {
		for (size_t j = i + 1; j < getSize(); j++) {
			if (vv_[i].isApprox(vv_[j], 0.0001)) { std::cerr << "Duplicate Yukawa solution!"; return true; }
		}
	}

	return false;

}

void Yukawa::splitPairsintoVectors(Yukawa& ys1, Yukawa& ys2) const {

	assert(solution_ == yo::PAIR);
	ys1.solution_ = ys2.solution_ = yo::VECTOR;

	bool con = false;
	ys1.reserve(getSize() / 2);
	ys2.reserve(getSize() / 2);
	for (auto i : vv_) {
		con = !con;
		if (con) { ys1.push_back(i); continue; }
		ys2.push_back(i);
	}

}

bool Yukawa::isUniqueVector(const MyVector& mv, const std::string& s) {

	assert(solution_ == yo::VECTOR);

	for (size_t i = 0; i < vv_.size(); i++) {
		if (vv_[i].isApprox(mv, 0.0001)) {
			if (vs_[i].find(s) == std::string::npos) vs_[i] += s;
			return false;
		}
	}

	vv_.push_back(mv); vs_.push_back(s);
	return true;

}

void findUniqueVectors(std::vector<Yukawa>& vys, Yukawa& ys, const std::string& s) {

	ys.solution_ = yo::VECTOR;

	for (auto i : vys) {
		for (auto j : i) ys.isUniqueVector(j, s);
	}

}

bool Yukawa::isUniquePair(const MyVector& mv1, const MyVector& mv2, const std::string& s) {

	assert(solution_ == yo::PAIR);

	for (size_t i = 0; i < vv_.size(); i += 2) {
		if (vv_[i].isApprox(mv1, 0.0001) && vv_[i + 1].isApprox(mv2, 0.0001)) {
			if (vs_[i / 2].find(s) == std::string::npos) vs_[i / 2] += s;
			return false;
		}
	}

	vv_.push_back(mv1); vv_.push_back(mv2); vs_.push_back(s);
	return true;

}

void findUniquePairs(std::vector<Yukawa>& vys, Yukawa& ys, const std::string& s) {

	ys.solution_ = yo::PAIR;

	for (size_t i = 0; i < vys.size() - 1; i += 2) {
		if ((vys[i].getSize() == 1) && (vys[i + 1].getSize() == 1)) { ys.isUniquePair(vys[i][0], vys[i + 1][0], s); continue; }
		for (size_t j = 0; j < vys[i].getSize(); j++) {
			for (size_t k = 0; k < vys[i + 1].getSize(); k++) { ys.isUniquePair(vys[i][j], vys[i + 1][k], s); }
		}
	}

}

bool Yukawa::isUniqueMatrix(const MyVector& mv, const std::string& s) {

	assert(solution_ == yo::VECTOR);

	std::vector<Eigen::Matrix3cd> vm = convertSolutionstoMassMatrices();
	Eigen::Matrix3cd m = mv.getMassMatrix() * mv.getMassMatrix().adjoint();

	for (size_t i = 0; i < vm.size(); i++) {
		if (vm[i].isApprox(m, 0.0001)) {
			vs_[i] += s + " ";
			return false;
		}
	}

	vv_.push_back(mv); vs_.push_back(s + " ");
	return true;

}

void findUniqueMatrices(const Yukawa& ysu_v, Yukawa& ysu_m) {

	ysu_m.solution_ = ysu_v.solution_;

	for (size_t i = 0; i < ysu_v.getSize(); i++) ysu_m.isUniqueMatrix(ysu_v[i], std::to_string(i + 1));

}

bool Yukawa::isUniqueMatrixPair(const MyVector& mv1, const MyVector& mv2, const std::string& s) {

	assert(solution_ == yo::PAIR);

	std::vector<Eigen::Matrix3cd> vm = convertSolutionstoMassMatrices();
	Eigen::Matrix3cd m1 = mv1.getMassMatrix() * mv1.getMassMatrix().adjoint();
	Eigen::Matrix3cd m2 = mv2.getMassMatrix() * mv2.getMassMatrix().adjoint();

	for (size_t i = 0; i < vm.size(); i += 2) {
		if (vm[i].isApprox(m1, 0.0001) && vm[i + 1].isApprox(m2, 0.0001)) {
			vs_[i / 2] += s + " ";
			return false;
		}
	}

	vv_.push_back(mv1); vv_.push_back(mv2); vs_.push_back(s + " ");
	return true;

}

void findUniqueMatrixPairs(const Yukawa& ysu_p, Yukawa& ysu_m) {

	ysu_m.solution_ = ysu_p.solution_;

	for (size_t i = 0; i < ysu_p.getSize(); i += 2) {
		ysu_m.isUniqueMatrixPair(ysu_p[i], ysu_p[i + 1], std::to_string((i / 2) + 1));
	}

}

void findMassRatio(const Yukawa& ys, double step, Yukawa& ys_m, std::vector<std::vector<double>>& vvd) {

	std::vector<double> vd;

	if ((int)ys.solution_ == 1) findUniqueMatrices(ys, ys_m);
	if ((int)ys.solution_ == 2) findUniqueMatrixPairs(ys, ys_m);
	vvd.reserve(ys_m.getSize());

	for (auto i : ys_m) {
		i.getMassRatio(step, vd);
		vvd.push_back(vd);
		vd.clear();
	}

}

void Yukawa::findSolutionsWithPhases(Yukawa& ys) const {

	ys.solution_ = solution_;

	if ((int)solution_ == 1) {
		for (size_t i = 0; i < getSize(); i++) ys.isUniqueVector(vv_[i].setAllNonZeroElementstoPhases());
	}

	if ((int)solution_ == 2) {
		for (size_t i = 0; i < getSize(); i += 2) ys.isUniquePair(vv_[i].setAllNonZeroElementstoPhases(), vv_[i + 1].setAllNonZeroElementstoPhases());
	}

}

void Yukawa::findSolutionsWithUnityElements(Yukawa& ys1, Yukawa& ys2) const {

	findSolutionsWithPhases(ys1);
	ys2.solution_ = solution_;

	if ((int)solution_ == 1) {
		std::vector<std::complex<double>> vcd;
		for (size_t i = 0; i < ys1.getSize(); i++) {
			(ys1[i].getPhases()).swap(vcd);
			ys2.isUniqueVector(ys1[i].setAllNonZeroElementsto1(), convertPhasestoString(vcd));
		}
	}

	if ((int)solution_ == 2) {
		std::vector<std::complex<double>> vcd1, vcd2;
		for (size_t i = 0; i < ys1.getSize(); i += 2) {
			(ys1[i].getPhases()).swap(vcd1); (ys1[i + 1].getPhases()).swap(vcd2);
			ys2.isUniquePair(ys1[i].setAllNonZeroElementsto1(), ys1[i + 1].setAllNonZeroElementsto1(), convertPhasestoString(vcd1, vcd2));
		}
	}

}

void Yukawa::printToFile(const std::string& s) {

	std::fstream ofile;
	int k;

	ofile = fileOpener("outputs/", s, std::ios::out);

	for (size_t i = 0; i < vv_.size(); i += (int)solution_) {
		k = 0;
		ofile << (i / (int)solution_) + 1 << std::endl;
		while (k < (int)solution_) {
			vv_[i + k].setFirstNonZeroElementto1(); vv_[i + k].setNumericZerotoActualZero();
			ofile << vv_[i + k] << std::endl;
			k++;
		}
		ofile << vs_[i / (int)solution_] << std::endl;
	}

}

void Yukawa::printToFile(const std::string& s, std::vector<std::vector<double>> vvd) {

	std::fstream ofile;
	std::vector<Eigen::Matrix3cd> vm;
	MyMatrix mm;

	ofile = fileOpener("outputs/", s, std::ios::out);
	(convertSolutionstoMassMatrices()).swap(vm);
	
	for (size_t i = 0; i < vm.size(); i++) {
		mm = MyMatrix(vm[i]); mm.setNumericZerotoActualZero();
		if ((int)solution_ == 1) ofile << i + 1 << " { " << vs_[i] << "} " << std::endl;
		if ((int)solution_ == 2 && !(i % (int)solution_)) ofile << (i / 2) + 1 << " { " << vs_[i / 2] << "} " << std::endl;
		ofile << mm << std::endl << std::endl;
		if (vvd.size()) printMassRatio(vvd[i], ofile);
	}

	ofile.close();

}

MyVector& Yukawa::operator[](size_t index) {

	if (index > getSize()) {
		std::cerr << "Vector subscript out of range!" << std::endl;
		exit(EXIT_FAILURE);
	}
	return vv_[index];

}

const MyVector& Yukawa::operator[](size_t index) const {

	if (index > getSize()) {
		std::cerr << "Vector subscript out of range!" << std::endl;
		exit(EXIT_FAILURE);
	}
	return vv_[index];

}

std::ostream& operator<<(std::ostream& os, const Yukawa& y) {

	for (auto i : y.vv_) os << i << std::endl;
	return os;

}