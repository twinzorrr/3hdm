#pragma once
#ifndef YUKAWA_H_
#define YUKAWA_H_

template<typename T> class MyVector;
template<typename T> class Group;
template<typename T> using Real = typename Eigen::NumTraits<T>::Real;
using yo = enum class yukawa_option { UNDEFINED, VECTOR, PAIR };


template<typename T> void findUniqueVectors(std::vector<Yukawa<T>>&, Yukawa<T>&, const std::string&);
template<typename T> void findUniquePairs(std::vector<Yukawa<T>>&, Yukawa<T>&, const std::string&);
template<typename T> void findUniqueMatrices(const Yukawa<T>&, Yukawa<T>&);
template<typename T> void findUniqueMatrixPairs(const Yukawa<T>&, Yukawa<T>&);
template<typename T> void findMassRatio(const Yukawa<T>&, Real<T>, Yukawa<T>&, std::vector<std::vector<Real<T>>>&);

template<typename _scalar>
class Yukawa
{
public:
	Yukawa();
	Yukawa(yo, const std::vector<MyVector<_scalar>>&, const std::vector<std::string>& vs = {});
	Yukawa(const Yukawa&) = default;
	Yukawa(Yukawa&&) = default;
	size_t getSize() const { return vv_.size(); }
	void reserve(size_t i) { vv_.reserve(i); }
	void emplace_back(const Eigen::Matrix<_scalar, Eigen::Dynamic, 1>& v) { vv_.emplace_back(v); }
	void push_back(const MyVector<_scalar>& mv) { vv_.push_back(mv); }
	void clear() { vv_.clear(); vs_.clear(); }
	void swap(Yukawa& y) { y.vv_.swap(vv_); y.vs_.swap(vs_); }
	void setNumericZerotoActualZero();
	void setFirstNonZeroElementto1();
	bool errorGuard(const Group<_scalar>&) const;
	void splitPairsintoVectors(Yukawa&, Yukawa&) const;
	bool isUniqueVector(const MyVector<_scalar>&, const std::string& s = "");
	friend void findUniqueVectors<_scalar>(std::vector<Yukawa>&, Yukawa&, const std::string&);
	bool isUniquePair(const MyVector<_scalar>&, const MyVector<_scalar>&, const std::string& s = "");
	friend void findUniquePairs<_scalar>(std::vector<Yukawa>&, Yukawa&, const std::string&);
	bool isUniqueMatrix(const MyVector<_scalar>&, const std::string& s = "");
	friend void findUniqueMatrices<_scalar>(const Yukawa&, Yukawa&);
	bool isUniqueMatrixPair(const MyVector<_scalar>&, const MyVector<_scalar>&, const std::string&);
	friend void findUniqueMatrixPairs<_scalar>(const Yukawa&, Yukawa&);
	friend void findMassRatio<_scalar>(const Yukawa&, Real<_scalar>, Yukawa&, std::vector<std::vector<Real<_scalar>>>&);
	void findSolutionsWithPhases(Yukawa&) const;
	void findSolutionsWithUnityElements(Yukawa&, Yukawa&) const;
	void printToFile(const std::string&);
	void printToFile(const std::string&, const std::vector<std::vector<Real<_scalar>>>& vvsc);
	auto begin() { return vv_.begin(); }
	auto end() { return vv_.end(); }
	auto begin() const { return vv_.begin(); }
	auto end() const { return vv_.end(); }
	MyVector<_scalar>& operator[](size_t);
	const MyVector<_scalar>& operator[](size_t) const;
	friend std::ostream& operator<<(std::ostream& os, const Yukawa& y) { for (auto i : y.vv_) os << i << std::endl; return os; }
	Yukawa& operator=(const Yukawa& other) { if (this != &other) { solution_ = other.solution_; vs_ = other.vs_; vv_ = other.vv_; } return *this; }
	Yukawa& operator=(Yukawa&& other) noexcept { if (this != &other) { solution_ = std::move(other.solution_); vs_ = std::move(other.vs_); vv_ = std::move(other.vv_); } return *this; }
private:
	std::vector<Eigen::Matrix<_scalar, 3, 3>> convertSolutionstoMassMatrices() const;
	void printMassRatio(const std::vector<Real<_scalar>>&, std::fstream&) const;
	yo solution_;
	std::vector<std::string> vs_;
	std::vector<MyVector<_scalar>> vv_;
};

template<typename _scalar>
std::vector<Eigen::Matrix<_scalar, 3, 3>> Yukawa<_scalar>::convertSolutionstoMassMatrices() const {

	Eigen::Matrix<_scalar, 3, 3> m;
	std::vector<Eigen::Matrix<_scalar, 3, 3>> vm;

	for (size_t i = 0; i < vv_.size(); i++) {
		m = vv_[i].getMassMatrix() * vv_[i].getMassMatrix().adjoint();
		vm.push_back(m);
	}

	return vm;

}

template<typename _scalar>
void Yukawa<_scalar>::printMassRatio(const std::vector<Real<_scalar>>& vsc, std::fstream& file) const {

	file << "Muon electron mass ratio" << std::endl << " min: " << vsc[0] << " max: " << vsc[1] << std::endl;
	file << "Tau electron mass ratio" << std::endl << " min: " << vsc[2] << " max: " << vsc[3] << std::endl;
	file << "Tau muon mass ratio" << std::endl << " min: " << vsc[4] << " max: " << vsc[5] << std::endl;
	file << std::endl << std::endl;

}

template<typename _scalar>
Yukawa<_scalar>::Yukawa() : solution_(yo::UNDEFINED) {}

template<typename _scalar>
Yukawa<_scalar>::Yukawa(yo solution, const std::vector<MyVector<_scalar>>& vmv, const std::vector<std::string>& vs)
	: solution_(solution), vs_(vs), vv_(vmv)
{
	for_each(vmv.begin(), vmv.end(), [](const MyVector<_scalar>& mv) { assert(mv.getNumberofElements() == 27); });
}

template<typename _scalar>
void Yukawa<_scalar>::setNumericZerotoActualZero() { for (auto& i : vv_) i.setNumericZerotoActualZero(); }

template<typename _scalar>
void Yukawa<_scalar>::setFirstNonZeroElementto1() { for (auto& i : vv_) i.setFirstNonZeroElementto1(); }

template<typename _scalar>
bool Yukawa<_scalar>::errorGuard(const Group<_scalar>& g) const {

	for (auto i : vv_) {
		if (!g.isEigenvector1(i)) { /*std::cerr << "Yukawa solution is not a eigenvector! ";*/ return true; }
		if (i.isApprox(MyVector<_scalar>(Eigen::Matrix<_scalar, 27, 1>::Zero(27)), (Real<_scalar>)0.01)) { /*std::cerr << "Yukawa solution is a zero vector! ";*/ return true; }
	}

	for (size_t i = 0; i < getSize() - 1; i++) {
		for (size_t j = i + 1; j < getSize(); j++) {
			if (vv_[i].isApprox(vv_[j], (Real<_scalar>)0.0001)) { /*std::cerr << "Duplicate Yukawa solution! ";*/ return true; }
		}
	}

	return false;

}

template<typename _scalar>
void Yukawa<_scalar>::splitPairsintoVectors(Yukawa<_scalar>& ys1, Yukawa<_scalar>& ys2) const {

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

template<typename _scalar>
bool Yukawa<_scalar>::isUniqueVector(const MyVector<_scalar>& mv, const std::string& s) {

	assert(solution_ == yo::VECTOR);

	for (size_t i = 0; i < vv_.size(); i++) {
		if (vv_[i].isApprox(mv, (Real<_scalar>)0.0001)) {
			if (vs_[i].find(s) == std::string::npos) vs_[i] += s;
			return false;
		}
	}

	vv_.push_back(mv); vs_.push_back(s);
	return true;

}

template<typename T>
void findUniqueVectors(std::vector<Yukawa<T>>& vys, Yukawa<T>& ys, const std::string& s) {

	ys.solution_ = yo::VECTOR;

	for (auto i : vys) {
		for (auto j : i) ys.isUniqueVector(j, s);
	}

}

template<typename _scalar>
bool Yukawa<_scalar>::isUniquePair(const MyVector<_scalar>& mv1, const MyVector<_scalar>& mv2, const std::string& s) {

	assert(solution_ == yo::PAIR);

	for (size_t i = 0; i < vv_.size(); i += 2) {
		if (vv_[i].isApprox(mv1, (Real<_scalar>)0.0001) && vv_[i + 1].isApprox(mv2, (Real<_scalar>)0.0001)) {
			if (vs_[i / 2].find(s) == std::string::npos) vs_[i / 2] += s;
			return false;
		}
	}

	vv_.push_back(mv1); vv_.push_back(mv2); vs_.push_back(s);
	return true;

}

template<typename T>
void findUniquePairs(std::vector<Yukawa<T>>& vys, Yukawa<T>& ys, const std::string& s) {

	ys.solution_ = yo::PAIR;

	for (size_t i = 0; i < vys.size() - 1; i += 2) {
		if ((vys[i].getSize() == 1) && (vys[i + 1].getSize() == 1)) { ys.isUniquePair(vys[i][0], vys[i + 1][0], s); continue; }
		for (size_t j = 0; j < vys[i].getSize(); j++) {
			for (size_t k = 0; k < vys[i + 1].getSize(); k++) { ys.isUniquePair(vys[i][j], vys[i + 1][k], s); }
		}
	}

}

template<typename _scalar>
bool Yukawa<_scalar>::isUniqueMatrix(const MyVector<_scalar>& mv, const std::string& s) {

	assert(solution_ == yo::VECTOR);

	std::vector<Eigen::Matrix<_scalar, 3, 3>> vm = convertSolutionstoMassMatrices();
	Eigen::Matrix<_scalar, 3, 3> m = mv.getMassMatrix() * mv.getMassMatrix().adjoint();

	for (size_t i = 0; i < vm.size(); i++) {
		if (vm[i].isApprox(m, (Real<_scalar>)0.0001)) {
			vs_[i] += s + " ";
			return false;
		}
	}

	vv_.push_back(mv); vs_.push_back(s + " ");
	return true;

}

template<typename T>
void findUniqueMatrices(const Yukawa<T>& ysu_v, Yukawa<T>& ysu_m) {

	ysu_m.solution_ = ysu_v.solution_;

	for (size_t i = 0; i < ysu_v.getSize(); i++) ysu_m.isUniqueMatrix(ysu_v[i], std::to_string(i + 1));

}

template<typename _scalar>
bool Yukawa<_scalar>::isUniqueMatrixPair(const MyVector<_scalar>& mv1, const MyVector<_scalar>& mv2, const std::string& s) {

	assert(solution_ == yo::PAIR);

	std::vector<Eigen::Matrix<_scalar, 3, 3>> vm = convertSolutionstoMassMatrices();
	Eigen::Matrix<_scalar, 3, 3> m1 = mv1.getMassMatrix() * mv1.getMassMatrix().adjoint();
	Eigen::Matrix<_scalar, 3, 3> m2 = mv2.getMassMatrix() * mv2.getMassMatrix().adjoint();

	for (size_t i = 0; i < vm.size(); i += 2) {
		if (vm[i].isApprox(m1, (Real<_scalar>)0.0001) && vm[i + 1].isApprox(m2, (Real<_scalar>)0.0001)) {
			vs_[i / 2] += s + " ";
			return false;
		}
	}

	vv_.push_back(mv1); vv_.push_back(mv2); vs_.push_back(s + " ");
	return true;

}

template<typename T>
void findUniqueMatrixPairs(const Yukawa<T>& ysu_p, Yukawa<T>& ysu_m) {

	ysu_m.solution_ = ysu_p.solution_;

	for (size_t i = 0; i < ysu_p.getSize(); i += 2) {
		ysu_m.isUniqueMatrixPair(ysu_p[i], ysu_p[i + 1], std::to_string((i / 2) + 1));
	}

}

template<typename T>
void findMassRatio(const Yukawa<T>& ys, Real<T> step, Yukawa<T>& ys_m, std::vector<std::vector<Real<T>>>& vvd) {

	std::vector<Real<T>> vd;

	if ((int)ys.solution_ == 1) { findUniqueMatrices(ys, ys_m); std::cout << std::endl << ys_m.getSize() << " unique mass matrices found!" << std::endl; }
	if ((int)ys.solution_ == 2) { findUniqueMatrixPairs(ys, ys_m); std::cout << std::endl << ys_m.getSize() / 2 << " unique pairs of mass matrices found!" << std::endl; }
	vvd.reserve(ys_m.getSize());

	for (auto i : ys_m) {
		i.getMassRatio(step, vd);
		vvd.push_back(vd);
		vd.clear();
	}

}

template<typename _scalar>
void Yukawa<_scalar>::findSolutionsWithPhases(Yukawa<_scalar>& ys) const {

	ys.solution_ = solution_;

	if ((int)solution_ == 1) {
		for (size_t i = 0; i < getSize(); i++) ys.isUniqueVector(vv_[i].setAllNonZeroElementstoPhases());
	}

	if ((int)solution_ == 2) {
		for (size_t i = 0; i < getSize(); i += 2) ys.isUniquePair(vv_[i].setAllNonZeroElementstoPhases(), vv_[i + 1].setAllNonZeroElementstoPhases());
	}

}

template<typename _scalar>
void Yukawa<_scalar>::findSolutionsWithUnityElements(Yukawa<_scalar>& ys1, Yukawa<_scalar>& ys2) const {

	findSolutionsWithPhases(ys1);
	ys2.solution_ = solution_;

	if ((int)solution_ == 1) {
		std::vector<_scalar> vsc;
		for (size_t i = 0; i < ys1.getSize(); i++) {
			(ys1[i].getPhases()).swap(vsc);
			ys2.isUniqueVector(ys1[i].setAllNonZeroElementsto1(), basic::convertPhasestoString(vsc));
		}
	}

	if ((int)solution_ == 2) {
		std::vector<_scalar> vsc1, vsc2;
		for (size_t i = 0; i < ys1.getSize(); i += 2) {
			(ys1[i].getPhases()).swap(vsc1); (ys1[i + 1].getPhases()).swap(vsc2);
			ys2.isUniquePair(ys1[i].setAllNonZeroElementsto1(), ys1[i + 1].setAllNonZeroElementsto1(), basic::convertPhasestoString(vsc1, vsc2));
		}
	}

}

template<typename _scalar>
void Yukawa<_scalar>::printToFile(const std::string& s) {

	std::fstream ofile;
	int k;

	ofile = basic::fileOpener("outputs/", s, std::ios::out);

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

template<typename _scalar>
void Yukawa<_scalar>::printToFile(const std::string& s, const std::vector<std::vector<Real<_scalar>>>& vvsc) {

	std::fstream ofile;
	std::vector<Eigen::Matrix<_scalar, 3, 3>> vm;
	MyMatrix<_scalar> mm;

	ofile = basic::fileOpener("outputs/", s, std::ios::out);
	(convertSolutionstoMassMatrices()).swap(vm);

	for (size_t i = 0; i < vm.size(); i++) {
		mm = MyMatrix<_scalar>(vm[i]); mm.setNumericZerotoActualZero();
		if ((int)solution_ == 1) ofile << i + 1 << " { " << vs_[i] << "} " << std::endl;
		if ((int)solution_ == 2 && !(i % (int)solution_)) ofile << (i / 2) + 1 << " { " << vs_[i / 2] << "} " << std::endl;
		ofile << mm << std::endl << std::endl;
		if (vvsc.size()) printMassRatio(vvsc[i], ofile);
	}

	ofile.close();

}

template<typename _scalar>
MyVector<_scalar>& Yukawa<_scalar>::operator[](size_t index) {

	if (index > getSize()) {
		std::cerr << "Vector subscript out of range!" << std::endl;
		exit(EXIT_FAILURE);
	}
	return vv_[index];

}

template<typename _scalar>
const MyVector<_scalar>& Yukawa<_scalar>::operator[](size_t index) const {

	if (index > getSize()) {
		std::cerr << "Vector subscript out of range!" << std::endl;
		exit(EXIT_FAILURE);
	}
	return vv_[index];

}

#endif // !YUKAWA_H_