#include <iostream>
#include <fstream>
#include "basicf.h"
#include "myvector.h"
#include "mymatrix.h"
#include "group.h"
#include "yukawa.h"


vector<Matrix3cd> Yukawa::convertSolutionstoMassMatrices() const {

	Matrix3cd m;
	vector<Matrix3cd> vm;

	for (size_t i = 0; i < vv_.size(); i++) {
		m = vv_[i].getMassMatrix() * vv_[i].getMassMatrix().adjoint();
		vm.push_back(m);
	}

	return vm;

}

void Yukawa::printMassRatio(vector<double>& vd, fstream& file) const {

	file << "Muon electron mass ratio" << endl << " min: " << vd[0] << " max: " << vd[1] << endl;
	file << "Tau electron mass ratio" << endl << " min: " << vd[2] << " max: " << vd[3] << endl;
	file << "Tau muon mass ratio" << endl << " min: " << vd[4] << " max: " << vd[5] << endl;
	file << endl << endl;

}

Yukawa::Yukawa() {}

Yukawa::Yukawa(vector<MyVector> vmv) { vv_ = vmv; }

void Yukawa::setNumericZerotoActualZero() { for (auto& i : vv_) i.setNumericZerotoActualZero(); }

void Yukawa::setFirstNonZeroElementto1() { for (auto& i : vv_) i.setFirstNonZeroElementto1(); }

bool Yukawa::errorGuard(const Group& g) const {

	for (auto i : vv_) {
		if (!g.isEigenvector1(i)) { cerr << "Yukawa solution is not a eigenvector!" << endl; return true; }
		if (i.isApprox(MyVector(VectorXcd::Zero(27)), 0.01)) { cerr << "Yukawa solution is a zero vector"; return true; }
	}

	for (size_t i = 0; i < getSize() - 1; i++) {
		for (size_t j = i + 1; j < getSize(); j++) {
			if (vv_[i].isApprox(vv_[j], 0.0001)) { cerr << "Duplicate Yukawa solution!"; return true; }
		}
	}

	return false;

}

void Yukawa::splitPairsintoVectors(Yukawa& ys1, Yukawa& ys2) const {

	bool con = false;
	ys1.reserve(getSize() / 2);
	ys2.reserve(getSize() / 2);
	for (auto i : vv_) {
		con = !con;
		if (con) { ys1.push_back(i); continue; }
		ys2.push_back(i);
	}

}

bool Yukawa::isUniqueVector(const MyVector& mv, const string& s) {

	for (size_t i = 0; i < vv_.size(); i++) {
		if (vv_[i].isApprox(mv, 0)) {
			if (vs_[i].find(s) == string::npos) vs_[i] += s;
			return false;
		}
	}

	vv_.push_back(mv); vs_.push_back(s);
	return true;

}

void findUniqueVectors(vector<Yukawa>& vys, Yukawa& ys, const string& s) {

	for (auto i : vys) {
		for (auto j : i) ys.isUniqueVector(j, s);
	}

}

bool Yukawa::isUniquePair(const MyVector& mv1, const MyVector& mv2, const string& s) {

	for (size_t i = 0; i < vv_.size(); i += 2) {
		if (vv_[i].isApprox(mv1, 0.0001) && vv_[i + 1].isApprox(mv2, 0.0001)) {
			if (vs_[i / 2].find(s) == string::npos) vs_[i / 2] += s;
			return false;
		}
	}

	vv_.push_back(mv1); vv_.push_back(mv2); vs_.push_back(s);
	return true;

}

void findUniquePairs(vector<Yukawa>& vys, Yukawa& ys, const string& s) {

	for (size_t i = 0; i < vys.size() - 1; i += 2) {
		if ((vys[i].getSize() == 1) && (vys[i + 1].getSize() == 1)) { ys.isUniquePair(vys[i][0], vys[i + 1][0], s); continue; }
		for (size_t j = 0; j < vys[i].getSize(); j++) {
			for (size_t k = 0; k < vys[i + 1].getSize(); k++) { ys.isUniquePair(vys[i][j], vys[i + 1][k], s); }
		}
	}

}

bool Yukawa::isUniqueMatrix(const MyVector& mv, const string& s) {

	vector<Matrix3cd> vm = convertSolutionstoMassMatrices();
	Matrix3cd m = mv.getMassMatrix() * mv.getMassMatrix().adjoint();

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

	for (size_t i = 0; i < ysu_v.getSize(); i++) ysu_m.isUniqueMatrix(ysu_v[i], to_string(i + 1));

}

void findMassRatio(const Yukawa& ysu_v, double step, Yukawa& ysu_m, vector<vector<double>>& vvd) {

	vector<double> vd;

	findUniqueMatrices(ysu_v, ysu_m);
	vvd.reserve(ysu_m.getSize());

	for (auto i : ysu_m) {
		i.getMassRatio(step, vd);
		vvd.push_back(vd);
		vd.clear();
	}

}

void Yukawa::findSolutionsWithPhases(Yukawa& ys) const {

	for (size_t i = 0; i < getSize(); i += 2) {
		ys.isUniquePair(vv_[i].setAllNonZeroElementstoPhases(), vv_[i + 1].setAllNonZeroElementstoPhases());
	}

}

void Yukawa::findSolutionsWithUnityElements(Yukawa& ys1, Yukawa& ys2) const {

	vector<complex<double>> vcd1, vcd2;

	findSolutionsWithPhases(ys1);

	for (size_t i = 0; i < ys1.getSize(); i += 2) {
		(ys1[i].getPhases()).swap(vcd1); (ys1[i + 1].getPhases()).swap(vcd2);
		ys2.isUniquePair(ys1[i].setAllNonZeroElementsto1(), ys1[i + 1].setAllNonZeroElementsto1(), convertPhasestoString(vcd1, vcd2));
	}

}

void Yukawa::printToFile(const string& s, po e) {

	fstream ofile;
	int k;

	ofile = fileOpener("outputs/", s, ios::out);

	for (size_t i = 0; i < vv_.size(); i += (int)e) {
		k = 0;
		ofile << (i / (int)e) + 1 << endl;
		while (k < (int)e) {
			vv_[i + k].setFirstNonZeroElementto1(); vv_[i + k].setNumericZerotoActualZero();
			ofile << vv_[i + k] << endl;
			k++;
		}
		ofile << vs_[i / (int)e] << endl;
	}

}

void Yukawa::printToFile(const string& s, vector<vector<double>> vvd) {

	fstream ofile;
	vector<Matrix3cd> vm;
	MyMatrix mm;

	ofile = fileOpener("outputs/", s, ios::out);
	(convertSolutionstoMassMatrices()).swap(vm);
	for (size_t i = 0; i < vm.size(); i++) {
		mm = MyMatrix(vm[i]); mm.setNumericZerotoActualZero();
		ofile << i + 1 << " { " << vs_[i] << "} " << endl << mm << endl << endl;
		if (vvd.size()) printMassRatio(vvd[i], ofile);
	}
	ofile.close();

}

MyVector& Yukawa::operator[](size_t index) {

	if (index > getSize()) {
		cerr << "Vector subscript out of range!" << endl;
		exit(EXIT_FAILURE);
	}
	return vv_[index];

}

const MyVector& Yukawa::operator[](size_t index) const {

	if (index > getSize()) {
		cerr << "Vector subscript out of range!" << endl;
		exit(EXIT_FAILURE);
	}
	return vv_[index];

}

ostream& operator<<(ostream& os, const Yukawa& y) {

	for (auto i : y.vv_) os << i << endl;
	return os;

}