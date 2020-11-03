#include <iostream>
#include <fstream>
#include "basicf.h"
#include "myvector.h"
#include "mymatrix.h"
#include "group.h"
#include "yukawa.h"


Group::Group() { group_ = "No Group"; nor_ = 0; nog_ = 0; }

Group::Group(vector<MyMatrix> vmm, string group, size_t nor) { vm_ = vmm; group_ = group; nor_ = nor; if (nor) nog_ = vmm.size() / nor; else nog_ = 0; }

void Group::setNumericZerotoActualZero() { for (auto& i : vm_) i.setNumericZerotoActualZero(); }

bool Group::isEigenvector1(const MyVector& mv) const {

	for (auto i : vm_) {
		if (!i.isEigenvector1(mv)) return false;
	}
	return true;

}

Group Group::findKroneckerProduct3(size_t a, size_t b, size_t c, const string& s) const {

	vector<MyMatrix> vmm;
	vmm.reserve(nog_);
	for (size_t i = 0; i < nog_; i++) { vmm.push_back(getKroneckerProduct3(vm_[a + i], vm_[b + i], vm_[c + i], s)); }
	return Group(vmm, group_, nor_);

}

Group Group::findEigenvectors1() const {

	MyMatrix mm;
	vector<MyMatrix> vmm;
	vmm.reserve(getSize());

	for (auto i : vm_) {
		mm = i.getEigenvectors1();
		if (mm.getNumberofCols() == 0) { vmm.clear(); break; }
		vmm.push_back(mm);
	}

	return Group(vmm, group_, nor_);

}

MyMatrix Group::findIntersectionBasis() const {

	size_t n;
	MyMatrix mm;

	n = 2;
	mm = getIntersectionBasis(vm_[0], vm_[1]);
	while (((mm.extractYukawaSolution()[0]).norm() != 0) && (n < vm_.size())) { mm = getIntersectionBasis(mm, vm_[n]); n++; }

	return mm;

}

vector<Yukawa> Group::findSolutions(const string& s, fstream& file, vector<string>* p_vs) const {

	string sout;
	Group kp, es;
	MyMatrix ibs;
	vector<Yukawa> vys;
	Yukawa ys;
	int nofs = 0;

	file << group_ << endl;
	for (size_t a = 0; a < getSize(); a += nog_) {
		for (size_t b = 0; b < getSize(); b += nog_) {
			for (size_t c = 0; c < getSize(); c += nog_) {
				sout = to_string((a / nog_) + 1) + "x" + to_string((b / nog_) + 1) + "x" + to_string((c / nog_) + 1);
				(findKroneckerProduct3(a, b, c, s)).swap(kp);
				(kp.findEigenvectors1()).swap(es);
				if (es.getSize() == 0) { flash(group_ + " " + sout + " No eigenvectors!", 40, 1); continue; }
				ibs = es.findIntersectionBasis();
				if ((ibs.extractYukawaSolution()[0]).norm() == 0) { flash(group_ + " " + sout + " No eigensubspace!", 40, 1); continue; }
				(ibs.extractYukawaSolution()).swap(ys); ys.setFirstNonZeroElementto1(); ys.setNumericZerotoActualZero();
				if (ys.errorGuard(kp)) { exit(EXIT_FAILURE); }
				flash(group_ + " " + sout + " Solution found!", 40, 1);
				vys.push_back(ys);
				p_vs->push_back(to_string(a) + " " + to_string(b));
				file << sout << " Yes" << endl << ys;
				nofs++;
			}
		}
	}
	cout << group_ << " Number of solutions: " << nofs << " (" + s + ") " << endl;
	file << endl;

	return vys;

}

vector<Yukawa> Group::findPairSolutions(vector<Yukawa>& vys1, vector<Yukawa>& vys2, fstream& file1, fstream& file2, fstream& file3) const {

	vector<Yukawa> vys;
	vector<string> vs1, vs2;
	bool justonce = 1;

	vys1 = findSolutions("charged", file1, &vs1);
	vys2 = findSolutions("dirac", file2, &vs2);

	for (size_t i = 0; i < vs1.size(); i++) {
		for (size_t j = 0; j < vs2.size(); j++) {
			if (vs1[i] == vs2[j]) {
				if (justonce) { file3 << group_ << endl; justonce = 0; }
				file3 << "{ ";
				for (auto k : vys1[i]) { k.setFirstNonZeroElementto1(); k.setNumericZerotoActualZero(); file3 << k << " "; }
				file3 << "} " << endl << "{ ";
				for (auto l : vys2[j]) { l.setFirstNonZeroElementto1(); l.setNumericZerotoActualZero(); file3 << l << " "; }
				file3 << "}" << endl << endl;
				vys.push_back(vys1[i]); vys.push_back(vys2[j]);
			}
		}
	}

	return vys;

}

MyMatrix& Group::operator[](size_t index) {

	if (index > getSize()) {
		cerr << "Vector subscript out of range!" << endl;
		exit(EXIT_FAILURE);
	}
	return vm_[index];

}

const MyMatrix& Group::operator[](size_t index) const {

	if (index > getSize()) {
		cerr << "Vector subscript out of range!" << endl;
		exit(EXIT_FAILURE);
	}
	return vm_[index];

}

ostream& operator<<(ostream& os, const Group& g) {

	for (auto i : g.vm_) os << i << endl;
	return os;

}