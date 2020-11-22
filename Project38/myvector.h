#pragma once
#ifndef MYVECTOR_H_
#define MYVECTOR_H_

template<typename T> using Real = typename Eigen::NumTraits<T>::Real;


template<typename T> int getIndexofFirstNonZeroElement(const Eigen::Matrix<T, Eigen::Dynamic, 1>&);
template<typename T> void setNumericZerotoActualZero_(Eigen::Matrix<T, Eigen::Dynamic, 1>&);
template<typename T> void setFirstNonZeroElementto1_(Eigen::Matrix<T, Eigen::Dynamic, 1>&);

template<typename _scalar>
class MyVector
{
public:
	MyVector();
	explicit MyVector(const Eigen::Matrix<_scalar, Eigen::Dynamic, 1>&);
	int getNumberofElements() const { return size_; }
	void setNumericZerotoActualZero();
	void setFirstNonZeroElementto1();
	Real<_scalar> norm() const { return v_.norm(); }
	bool isApprox(const MyVector&, Real<_scalar>) const;
	Eigen::Matrix<_scalar, 3, 3> getMassMatrix(Real<_scalar> v1 = 1, Real<_scalar> v2 = 1, Real<_scalar> v3 = 1, Real<_scalar> xi2 = acos(-1) / 2, Real<_scalar> xi3 = acos(-1) / 2) const;
	void getMassRatio(Real<_scalar>, std::vector<Real<_scalar>>&) const;
	std::vector<_scalar> getPhases();
	MyVector setAllNonZeroElementstoPhases() const;
	MyVector setAllNonZeroElementsto1() const;
	friend MyVector operator*(const Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>& m, const MyVector& mv) { return MyVector<_scalar>(m * mv.v_); }
	friend std::ostream& operator<<(std::ostream& os, const MyVector& mv) { return os << mv.v_.format(fmtV); }
private:
	friend int getIndexofFirstNonZeroElement(const Eigen::Matrix<_scalar, Eigen::Dynamic, 1>&);
	friend void setNumericZerotoActualZero_(Eigen::Matrix<_scalar, Eigen::Dynamic, 1>&);
	friend void setFirstNonZeroElementto1_(Eigen::Matrix<_scalar, Eigen::Dynamic, 1>&);
	std::vector<Eigen::Matrix<_scalar, 3, 3>> getYukawaMatrices_() const;
	int size_;
	Eigen::Matrix<_scalar, Eigen::Dynamic, 1> v_;
};

template<typename T>
int getIndexofFirstNonZeroElement(const Eigen::Matrix<T, Eigen::Dynamic, 1>& v) {

	for (int i = 0; i < v.size(); i++) { if (!(abs(v(i)) < Real<T>(LDBL_EPSILON * 10.0L))) return i; }
	return 0;

}

template<typename T>
void setNumericZerotoActualZero_(Eigen::Matrix<T, Eigen::Dynamic, 1>& v) {

	for (int i = 0; i < v.size(); i++) {
		if (abs(real(v(i))) < (Real<T>)0.0001) v(i).real((Real<T>)0);
		if (abs(imag(v(i))) < (Real<T>)0.0001) v(i).imag((Real<T>)0);
	}

}

template<typename T>
void setFirstNonZeroElementto1_(Eigen::Matrix<T, Eigen::Dynamic, 1>& v) {

	for (int i = 0; i < v.size(); i++) {
		if (!(abs(v(i)) < Real<T>(LDBL_EPSILON * 10.0L))) { v = ((Real<T>)1 / v(i)) * v; break; }
	}

}

template<typename _scalar>
std::vector<Eigen::Matrix<_scalar, 3, 3>> MyVector<_scalar>::getYukawaMatrices_() const {

	Eigen::Matrix<_scalar, 3, 3> h;
	std::vector<Eigen::Matrix<_scalar, 3, 3>> vh;
	vh.reserve(3);

	for (int i = 0; i < getNumberofElements(); i += 9) {
		h = Eigen::Map<const Eigen::Matrix<_scalar, 3, 3>>(v_.segment(i, 9).data());
		vh.push_back(h.transpose());
	}

	return vh;

}

template<typename _scalar>
MyVector<_scalar>::MyVector() : size_(0), v_(Eigen::Matrix<_scalar, 0, 1>::Zero(size_)) { assert(Eigen::NumTraits<_scalar>::IsComplex); }

template<typename _scalar>
MyVector<_scalar>::MyVector(const Eigen::Matrix<_scalar, Eigen::Dynamic, 1>& v) : size_(v.size()), v_(v) { assert(Eigen::NumTraits<_scalar>::IsComplex); }

template<typename _scalar>
void MyVector<_scalar>::setNumericZerotoActualZero() { setNumericZerotoActualZero_(v_); }

template<typename _scalar>
void MyVector<_scalar>::setFirstNonZeroElementto1() { setFirstNonZeroElementto1_(v_); }

template<typename _scalar>
bool MyVector<_scalar>::isApprox(const MyVector& mv, Real<_scalar> real) const { return v_.isApprox(mv.v_, real); }

template<typename _scalar>
Eigen::Matrix<_scalar, 3, 3> MyVector<_scalar>::getMassMatrix(Real<_scalar> v1, Real<_scalar> v2, Real<_scalar> v3, Real<_scalar> xi2, Real<_scalar> xi3) const {

	assert(v_.size() == 27);
	consts<Real<_scalar>> cs;
	return ((Real<_scalar>)(-1) / sqrt(2)) * (v1 * getYukawaMatrices_()[0] + exp(-cs.i * xi2) * v2 * getYukawaMatrices_()[1] + exp(-cs.i * xi3) * v3 * getYukawaMatrices_()[2]);

}

template<typename _scalar>
void MyVector<_scalar>::getMassRatio(Real<_scalar> step, std::vector<Real<_scalar>>& vn) const {

	consts<Real<_scalar>> cs;
	Real<_scalar> v1;
	Eigen::Matrix<_scalar, 3, 3> m;
	Eigen::Matrix<Real<_scalar>, 3, 3> diag;
	std::vector<Real<_scalar>> mume, mtme, mtmu;
	std::vector<std::vector<Real<_scalar>>> vvn;
	vn.reserve(6);

	for (Real<_scalar> v2 = 1; v2 < (Real<_scalar>)246; v2 += step) {
		for (Real<_scalar> v3 = v2; v3 < sqrt(pow((Real<_scalar>)246, (Real<_scalar>)2) - pow(v2, (Real<_scalar>)2)); v3 += step) {
			if ((v2 == (Real<_scalar>)0) & (v3 == (Real<_scalar>)0)) continue;
			for (Real<_scalar> xi2 = -1 * cs.pi; xi2 < cs.pi; xi2 += cs.pi / 2) {
				for (Real<_scalar> xi3 = -1 * cs.pi; xi3 < cs.pi; xi3 += cs.pi / 2) {
					v1 = sqrt(pow((Real<_scalar>)246, (Real<_scalar>)2) - pow(v2, (Real<_scalar>)2) - pow(v3, (Real<_scalar>)2));
					basic::flash(std::to_string(v1) + " " + std::to_string(v2) + " " + std::to_string(v3), 45, 1);
					m = getMassMatrix(v1, v2, v3, xi2, xi3); m = m * m.adjoint();
					Eigen::JacobiSVD<Eigen::Matrix<_scalar, 3, 3>> svd(m, Eigen::ComputeFullU);
					diag = svd.singularValues().asDiagonal();
					mume.emplace_back(sqrt(diag(1, 1) / diag(2, 2)));
					mtme.emplace_back(sqrt(diag(0, 0) / diag(2, 2)));
					mtmu.emplace_back(sqrt(diag(0, 0) / diag(1, 1)));
				}
			}
		}
	}
	std::cout << std::endl;

	vvn = { mume,mtme,mtmu };
	for (auto i : vvn) {
		auto minmax = minmax_element(begin(i), end(i));
		vn.insert(vn.end(), { *minmax.first,*minmax.second });
	}

}

template<typename _scalar>
std::vector<_scalar> MyVector<_scalar>::getPhases() {

	std::vector<_scalar> vsc;
	bool justonce = true;

	for (int i = 0; i < getNumberofElements(); i++) {
		if (!(abs(v_(i)) < (Real<_scalar>)LDBL_EPSILON * 10L)) {
			if (justonce) { justonce = false; continue; }
			vsc.push_back(v_(i));
		}
	}

	return vsc;

}

template<typename _scalar>
MyVector<_scalar> MyVector<_scalar>::setAllNonZeroElementstoPhases() const {

	assert(v_.size() % 9 == 0 && v_.size() > 9);

	std::vector<Eigen::Matrix<_scalar, Eigen::Dynamic, 1>> vv;
	Eigen::Matrix<_scalar, Eigen::Dynamic, 1> v;

	for (int i = 0; i < v_.size(); i += 9) vv.push_back(v_.segment(i, 9));
	std::sort(vv.begin(), vv.end(), [](Eigen::Matrix<_scalar, Eigen::Dynamic, 1> x, Eigen::Matrix<_scalar, Eigen::Dynamic, 1> y) { return getIndexofFirstNonZeroElement(x) < getIndexofFirstNonZeroElement(y); });
	for (size_t i = 0; i < vv.size(); i++) { v.conservativeResize((i + 1) * 9); v.tail(9) = vv[i]; }
	setFirstNonZeroElementto1_(v);
	return MyVector(v);

}

template<typename _scalar>
MyVector<_scalar> MyVector<_scalar>::setAllNonZeroElementsto1() const {

	Eigen::Matrix<_scalar, Eigen::Dynamic, 1> v = Eigen::Matrix<_scalar, Eigen::Dynamic, 1>::Zero(27);
	
	for (int i = 0; i < getNumberofElements(); i++) {
		if (!(abs(v_(i)) < (Real<_scalar>)LDBL_EPSILON * 10L)) {
			v(i).real((Real<_scalar>)1); v(i).imag((Real<_scalar>)0);
		}
	}
	
	return MyVector(v);

}

#endif // !MYVECTOR_H_