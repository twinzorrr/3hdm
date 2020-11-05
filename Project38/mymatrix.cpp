#include "Eigen/KroneckerProduct"
#include "myvector.h"
#include "mymatrix.h"
#include "yukawa.h"

static Eigen::IOFormat fmtM(6, 0, " ", "\n", "[", "]");


void setNumericZerotoActualZero_(Eigen::MatrixXcd& m) {

	for (int row = 0; row < m.rows(); row++) {
		for (int col = 0; col < m.cols(); col++) {
			if (abs(real(m(row, col))) < 0.0001) { m(row, col).real(0.0); }
			if (abs(imag(m(row, col))) < 0.0001) { m(row, col).imag(0.0); }
		}
	}

}

MyMatrix::MyMatrix() { row_ = 3; col_ = 3; m_ = Eigen::MatrixXcd::Zero(row_, col_); }

MyMatrix::MyMatrix(Eigen::MatrixXcd m) { row_ = m.rows(); col_ = m.cols(); m_ = m; }

void MyMatrix::setNumericZerotoActualZero() { setNumericZerotoActualZero_(m_); }

bool MyMatrix::isEigenvector1(const MyVector& mv) const {

	if ((m_ * mv).isApprox(mv, 0.0001)) return true;
	return false;

}

MyMatrix getKroneckerProduct3(const MyMatrix& mm1, const MyMatrix& mm2, const MyMatrix& mm3, const std::string& s) {

	Eigen::Matrix3cd H1, T1, H2, T3;
	Eigen::MatrixXcd KP, KPP;

	H1 = (mm1.m_).adjoint();
	T1 = (mm1.m_).transpose();
	H2 = (mm2.m_).adjoint();
	T3 = (mm3.m_).transpose();

	if (s == "charged") { KP = kroneckerProduct(H1, H2); KPP = kroneckerProduct(KP, T3); }
	if (s == "dirac") { KP = kroneckerProduct(T1, H2); KPP = kroneckerProduct(KP, T3); }

	return MyMatrix(KPP);

}

MyMatrix MyMatrix::getEigenvectors1() const {

	Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces;
	Eigen::MatrixXcd m;
	size_t j, n;
	bool con;

	ces.compute(m_);
	m.resize(27, (ces.eigenvectors()).size());
	j = 0; n = 0;

	for (int i = 0; i < (ces.eigenvalues()).size(); i++) {
		con = (abs(real(ces.eigenvalues()[i]) - 1.0) < 0.0001) & (abs(imag(ces.eigenvalues()[i])) < 0.0001);
		if (con) { m.col(j) = (ces.eigenvectors()).col(i); j++; n++; }
	}
	m.conservativeResize(27, n);

	return MyMatrix(m);

}

MyMatrix getIntersectionBasis(const MyMatrix& mm1, const MyMatrix& mm2) {

	Eigen::MatrixXcd m, mn, mm;
	Eigen::VectorXcd v;

	m.resize(27, mm1.col_ + mm2.col_);
	m << mm1.m_, -1.0 * mm2.m_;
	Eigen::FullPivLU<Eigen::MatrixXcd> lu(m);
	mn = lu.kernel();

	if (!lu.dimensionOfKernel()) { return MyMatrix(mn); }
	else {
		mm.resize(27, mn.cols());
		for (int i = 0; i < mn.cols(); i++) {
			v.setZero(27);
			for (int j = 0; j < mm1.col_; j++) {
				v += mn.col(i)(j) * (mm1.m_).col(j);
			}
			mm.col(i) = v;
		}
		return MyMatrix(mm);
	}

}

Yukawa MyMatrix::extractYukawaSolution() const {

	Yukawa ys;
	ys.reserve(col_);
	for (int col = 0; col < col_; col++) { ys.emplace_back(m_.col(col)); }
	return ys;

}

std::ostream& operator<<(std::ostream& os, const MyMatrix& mm) {

	return os << (mm.m_).format(fmtM);

}