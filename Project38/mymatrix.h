#pragma once
#ifndef MYMATRIX_H_
#define MYMATRIX_H_

template<typename T> class MyVector;
template<typename T> class MyMatrix;
template<typename T> class Yukawa;
template<typename T> using Real = typename Eigen::NumTraits<T>::Real;


template<typename T> void setNumericZerotoActualZero_(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>&);
template<typename T, typename U> MyMatrix<T> mycast(const MyMatrix<U>&, bool, bool);
template<typename T> MyMatrix<T> getKroneckerProduct3(const MyMatrix<T>&, const MyMatrix<T>&, const MyMatrix<T>&, const std::string&);
template<typename T> MyMatrix<T> getIntersectionBasis(const MyMatrix<T>&, const MyMatrix<T>&);

template<typename _scalar>
class MyMatrix
{
public:
	MyMatrix();
	explicit MyMatrix(const Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>&);
	MyMatrix(const MyMatrix&) = default;
	MyMatrix(MyMatrix&&) = default;
	int getNumberofRows() const { return row_; }
	int getNumberofCols() const { return col_; }
	void setNumericZerotoActualZero();
	bool isEigenvector1(const MyVector<_scalar>&) const;
	friend MyMatrix getKroneckerProduct3<_scalar>(const MyMatrix&, const MyMatrix&, const MyMatrix&, const std::string&);
	MyMatrix getEigenvectors1() const;
	friend MyMatrix getIntersectionBasis<_scalar>(const MyMatrix&, const MyMatrix&);
	Yukawa<_scalar> extractYukawaSolution() const;
	friend MyMatrix mycast(const MyMatrix&, bool, bool);
	const _scalar operator()(int row, int col) const { return m_(row, col); }
	friend std::ostream& operator<<(std::ostream& os, const MyMatrix& mm) { return os << (mm.m_).format(fmtM); }
	MyMatrix& operator=(const MyMatrix& other) { if (this != &other) { row_ = other.row_; col_ = other.col_; m_ = other.m_; } return *this; }
	MyMatrix& operator=(MyMatrix&& other) noexcept { if (this != &other) { row_ = std::move(other.row_); col_ = std::move(other.col_); m_ = std::move(other.m_); } return *this; }
private:
	friend void setNumericZerotoActualZero_(Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>&);
	int row_;
	int col_;
	Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic> m_;
};

template<typename T>
void setNumericZerotoActualZero_(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& m) {

	for (int row = 0; row < m.rows(); row++) {
		for (int col = 0; col < m.cols(); col++) {
			if (abs(real(m(row, col))) < (Real<T>)0.000001) { m(row, col).real((Real<T>)0); }
			if (abs(imag(m(row, col))) < (Real<T>)0.000001) { m(row, col).imag((Real<T>)0); }
		}
	}

}

template<typename _scalar>
MyMatrix<_scalar>::MyMatrix() : row_(0), col_(0), m_(Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>::Zero(row_, col_)) { assert(Eigen::NumTraits<_scalar>::IsComplex); }

template<typename _scalar>
MyMatrix<_scalar>::MyMatrix(const Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>& m) : row_(m.rows()), col_(m.cols()), m_(m) { assert(Eigen::NumTraits<_scalar>::IsComplex); }

template<typename _scalar>
void MyMatrix<_scalar>::setNumericZerotoActualZero() { setNumericZerotoActualZero_(m_); }

template<typename _scalar>
bool MyMatrix<_scalar>::isEigenvector1(const MyVector<_scalar>& mv) const {

	if ((m_ * mv).isApprox(mv, (Real<_scalar>)0.0001)) return true;
	return false;

}

template<typename T>
MyMatrix<T> getKroneckerProduct3(const MyMatrix<T>& mm1, const MyMatrix<T>& mm2, const MyMatrix<T>& mm3, const std::string& s) {

	Eigen::Matrix<T, 3, 3> H1, T1, H2, T3;
	Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> KP, KPP;

	H1 = (mm1.m_).adjoint();
	T1 = (mm1.m_).transpose();
	H2 = (mm2.m_).adjoint();
	T3 = (mm3.m_).transpose();

	if (s == "charged") { KP = kroneckerProduct(H1, H2); KPP = kroneckerProduct(KP, T3); }
	if (s == "dirac") { KP = kroneckerProduct(T1, H2); KPP = kroneckerProduct(KP, T3); }

	return MyMatrix<T>(KPP);

}

template<typename _scalar>
MyMatrix<_scalar> MyMatrix<_scalar>::getEigenvectors1() const {

	Eigen::ComplexEigenSolver<Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic>> ces;
	Eigen::Matrix<_scalar, Eigen::Dynamic, Eigen::Dynamic> m;
	size_t j, n;
	bool con;

	ces.compute(m_);
	m.resize(27, (ces.eigenvectors()).size());
	j = 0; n = 0;

	for (int i = 0; i < (ces.eigenvalues()).size(); i++) {
		con = (abs(real(ces.eigenvalues()[i]) - (Real<_scalar>)1) < (Real<_scalar>)0.0001) && (abs(imag(ces.eigenvalues()[i])) < (Real<_scalar>)0.0001);
		if (con) { m.col(j) = (ces.eigenvectors()).col(i); j++; n++; }
	}
	m.conservativeResize(27, n);

	return MyMatrix(m);

}

template<typename T>
MyMatrix<T> getIntersectionBasis(const MyMatrix<T>& mm1, const MyMatrix<T>& mm2) {

	Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> m, mn, mm;
	Eigen::Matrix<T, Eigen::Dynamic, 1> v;

	m.resize(27, mm1.col_ + mm2.col_);
	m << mm1.m_, (Real<T>)-1 * mm2.m_;
	Eigen::FullPivLU<Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>> lu(m);
	mn = lu.kernel();

	if (!lu.dimensionOfKernel()) { return MyMatrix<T>(mn); }
	else {
		mm.resize(27, mn.cols());
		for (int i = 0; i < mn.cols(); i++) {
			v.setZero(27);
			for (int j = 0; j < mm1.col_; j++) {
				v += mn.col(i)(j) * (mm1.m_).col(j);
			}
			mm.col(i) = v;
		}
		return MyMatrix<T>(mm);
	}

}

template<typename _scalar>
Yukawa<_scalar> MyMatrix<_scalar>::extractYukawaSolution() const {

	Yukawa<_scalar> ys;
	ys.reserve(col_);
	for (int col = 0; col < col_; col++) { ys.emplace_back(m_.col(col)); }
	return ys;

}

template<typename T, typename U>
MyMatrix<T> mycast(const MyMatrix<U>& mm, bool iscast1, bool iscast2) {

	Eigen::Matrix<T, 27, 27> nn;
	size_t l1, l2;
	std::string sreal, simag, s;
	std::stringstream ss;

	for (int i = 0; i < mm.getNumberofRows(); i++) {
		for (int j = 0; j < mm.getNumberofCols(); j++) {
			sreal = boost::lexical_cast<std::string>(mm(i, j).real());
			simag = boost::lexical_cast<std::string>(mm(i, j).imag());
			l1 = 10; l2 = 10;
			if (sreal[0] == '-') l1 = 11; if (simag[0] == '-') l2 = 11;
			std::string s1(sreal, 0, l1); std::string s2(simag, 0, l2);
			if (iscast1) { if (s2 == "0.64278760") s2 = "0.64278761"; if (s2 == "-0.64278760") s2 = "-0.64278761"; }
			if (iscast2) { if (s1 == "0.17364817") s1 = "0.17364818"; if (s1 == "-0.17364817") s1 = "-0.17364818"; }
			s = " (" + s1 + "," + s2 + ")";
			ss << s;
		}
	}

	load::detail::loadMatrix(ss, nn, 27);

	return MyMatrix<T>(nn);

}

#endif // !MYMATRIX_H_