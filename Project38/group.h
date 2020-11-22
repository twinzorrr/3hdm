#pragma once
#ifndef GROUP_H_
#define GROUP_H_

template<typename T> class MyVector;
template<typename T> class MyMatrix;
template<typename T> class Yukawa;


template<typename _scalar>
class Group
{
public:
	Group();
	explicit Group(const std::vector<MyMatrix<_scalar>>&, const std::string& group = "No Group", size_t nor = 0);
	size_t getSize() const { return vm_.size(); }
	void reserve(size_t i) { vm_.reserve(i); }
	void emplace_back(const Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>& m) { vm_.emplace_back(m); }
	void push_back(const MyMatrix<_scalar>& mm) { vm_.push_back(mm); }
	void clear() { vm_.clear(); }
	void swap(Group& g) { (g.vm_).swap(vm_); (g.group_).swap(group_); size_t tmp; tmp = g.nor_; g.nor_ = nor_; nor_ = tmp; }
	void setNumericZerotoActualZero();
	bool isEigenvector1(const MyVector<_scalar>&) const;
	Group findKroneckerProduct3(size_t, size_t, size_t, const std::string&) const;
	Group findEigenvectors1() const;
	MyMatrix<_scalar> findIntersectionBasis() const;
	std::vector<Yukawa<_scalar>> findSolutions(const std::string&, std::fstream&, std::vector<std::string>* p_vs = nullptr) const;
	std::vector<Yukawa<_scalar>> findPairSolutions(std::vector<Yukawa<_scalar>>&, std::vector<Yukawa<_scalar>>&, std::fstream&, std::fstream&, std::fstream&) const;
	auto begin() { return vm_.begin(); }
	auto end() { return vm_.end(); }
	auto begin() const { return vm_.begin(); }
	auto end() const { return vm_.end(); }
	MyMatrix<_scalar>& operator[](size_t);
	const MyMatrix<_scalar>& operator[](size_t) const;
	friend std::ostream& operator<<(std::ostream& os, const Group& g) { for (auto i : g.vm_) os << i << std::endl; return os; }
private:
	size_t nor_;
	size_t nog_;
	std::string group_;
	std::vector<MyMatrix<_scalar>> vm_;
};

template<typename _scalar>
Group<_scalar>::Group() : nor_(0), nog_(0), group_("No group") {}

template<typename _scalar>
Group<_scalar>::Group(const std::vector<MyMatrix<_scalar>>& vmm, const std::string& group, size_t nor)
	: nor_(nor), nog_((nor != 0) ? vmm.size() / nor : 0), group_(group), vm_(vmm)
{}

template<typename _scalar>
void Group<_scalar>::setNumericZerotoActualZero() { for (auto& i : vm_) i.setNumericZerotoActualZero(); }

template<typename _scalar>
bool Group<_scalar>::isEigenvector1(const MyVector<_scalar>& mv) const {

	for (auto i : vm_) {
		if (!i.isEigenvector1(mv)) return false;
	}
	return true;

}

template<typename _scalar>
Group<_scalar> Group<_scalar>::findKroneckerProduct3(size_t a, size_t b, size_t c, const std::string& s) const {

	std::vector<MyMatrix<_scalar>> vmm;
	vmm.reserve(nog_);
	for (size_t i = 0; i < nog_; i++) { vmm.push_back(getKroneckerProduct3(vm_[a + i], vm_[b + i], vm_[c + i], s)); }
	return Group<_scalar>(vmm, group_, nor_);

}

template<typename _scalar>
Group<_scalar> Group<_scalar>::findEigenvectors1() const {

	MyMatrix<_scalar> mm;
	std::vector<MyMatrix<_scalar>> vmm;
	vmm.reserve(getSize());

	for (auto i : vm_) {
		mm = i.getEigenvectors1();
		if (mm.getNumberofCols() == 0) { vmm.clear(); break; }
		vmm.push_back(mm);
	}

	return Group<_scalar>(vmm, group_, nor_);

}

template<typename _scalar>
MyMatrix<_scalar> Group<_scalar>::findIntersectionBasis() const {

	size_t n;
	MyMatrix<_scalar> mm;

	n = 2;
	mm = getIntersectionBasis(vm_[0], vm_[1]);
	while (((mm.extractYukawaSolution()[0]).norm() != (Real<_scalar>)0) && (n < vm_.size())) { mm = getIntersectionBasis(mm, vm_[n]); n++; }

	return mm;

}

template<typename _scalar>
std::vector<Yukawa<_scalar>> Group<_scalar>::findSolutions(const std::string& s, std::fstream& file, std::vector<std::string>* p_vs) const {

	std::string sout;
	Group<_scalar> kp, es;
	MyMatrix<_scalar> ibs;
	std::vector<Yukawa<_scalar>> vys;
	Yukawa<_scalar> ys;
	int nofs = 0;

	file << group_ << std::endl;
	for (size_t a = 0; a < getSize(); a += nog_) {
		for (size_t b = 0; b < getSize(); b += nog_) {
			for (size_t c = 0; c < getSize(); c += nog_) {
				sout = std::to_string((a / nog_) + 1) + "x" + std::to_string((b / nog_) + 1) + "x" + std::to_string((c / nog_) + 1);
				(findKroneckerProduct3(a, b, c, s)).swap(kp);
				(kp.findEigenvectors1()).swap(es);
				if (es.getSize() == 0) { basic::flash(group_ + " " + sout + " No eigenvectors!", 40, 1); continue; }
				ibs = es.findIntersectionBasis();
				if ((ibs.extractYukawaSolution()[0]).norm() == (Real<_scalar>)0) { basic::flash(group_ + " " + sout + " No eigensubspace!", 40, 1); continue; }
				(ibs.extractYukawaSolution()).swap(ys); ys.setFirstNonZeroElementto1(); ys.setNumericZerotoActualZero();
				if (ys.errorGuard(kp)) { exit(EXIT_FAILURE); }
				basic::flash(group_ + " " + sout + " Solution found!", 40, 1);
				vys.push_back(ys);
				if (p_vs) p_vs->push_back(std::to_string(a) + " " + std::to_string(b));
				file << sout << " Yes" << std::endl << ys;
				nofs++;
			}
		}
	}
	std::cout << group_ << " Number of solutions: " << nofs << " (" + s + ") " << std::endl;
	file << std::endl;

	return vys;

}

template<typename _scalar>
std::vector<Yukawa<_scalar>> Group<_scalar>::findPairSolutions(std::vector<Yukawa<_scalar>>& vys1, std::vector<Yukawa<_scalar>>& vys2, std::fstream& file1, std::fstream& file2, std::fstream& file3) const {

	std::vector<Yukawa<_scalar>> vys;
	std::vector<std::string> vs1, vs2;
	bool justonce = 1;

	vys1 = findSolutions("charged", file1, &vs1);
	vys2 = findSolutions("dirac", file2, &vs2);

	for (size_t i = 0; i < vs1.size(); i++) {
		for (size_t j = 0; j < vs2.size(); j++) {
			if (vs1[i] == vs2[j]) {
				if (justonce) { file3 << group_ << std::endl; justonce = 0; }
				file3 << "{ ";
				for (auto k : vys1[i]) { k.setFirstNonZeroElementto1(); k.setNumericZerotoActualZero(); file3 << k << " "; }
				file3 << "} " << std::endl << "{ ";
				for (auto l : vys2[j]) { l.setFirstNonZeroElementto1(); l.setNumericZerotoActualZero(); file3 << l << " "; }
				file3 << "}" << std::endl << std::endl;
				vys.push_back(vys1[i]); vys.push_back(vys2[j]);
			}
		}
	}

	return vys;

}

template<typename _scalar>
MyMatrix<_scalar>& Group<_scalar>::operator[](size_t index) {

	if (index > getSize()) {
		std::cerr << "Vector subscript out of range!" << std::endl;
		exit(EXIT_FAILURE);
	}
	return vm_[index];

}

template<typename _scalar>
const MyMatrix<_scalar>& Group<_scalar>::operator[](size_t index) const {

	if (index > getSize()) {
		std::cerr << "Vector subscript out of range!" << std::endl;
		exit(EXIT_FAILURE);
	}
	return vm_[index];

}

#endif // !GROUP_H_
