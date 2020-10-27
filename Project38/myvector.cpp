#include <iostream>
#include "basicf.h"
#include "myvector.h"

static IOFormat fmtV(6, DontAlignCols, " ", " ", "", "", "[", "]");
static const complex<double> i(0.0, 1.0);
static const double pi = acos(-1.0);


void setNumericZerotoActualZero_(VectorXcd& v) {

	for (int i = 0; i < v.size(); i++) {
		if (abs(real(v(i))) < 0.0001) v(i).real(0.0);
		if (abs(imag(v(i))) < 0.0001) v(i).imag(0.0);
	}

}

void setFirstNonZeroElementto1_(VectorXcd& v) {

	for (int i = 0; i < v.size(); i++) {
		if (!(abs(v(i)) < DBL_EPSILON * 10)) { v = (1.0 / v(i)) * v; break; }
	}

}

vector<Matrix3cd> MyVector::getYukawaMatrices_() const {

	Matrix3cd h;
	vector<Matrix3cd> vh;
	vh.reserve(3);

	for (int i = 0; i < getNumberofElements(); i += 9) {
		h = Map<Matrix<complex<double>, 3, 3>>((segment(i, 9)).array()).transpose();
		vh.push_back(h.transpose());
	}

	return vh;

}

MyVector::MyVector() { size_ = 3; v_ = VectorXcd::Zero(size_); }

MyVector::MyVector(VectorXcd v) { v_ = v; size_ = v.size(); }

void MyVector::setNumericZerotoActualZero() { setNumericZerotoActualZero_(v_); }

void MyVector::setFirstNonZeroElementto1() { setFirstNonZeroElementto1_(v_); }

bool MyVector::isApprox(const MyVector& mv, double d) const { return v_.isApprox(mv.v_, d); }

Matrix3cd MyVector::getMassMatrix(double v1, double v2, double v3, double xi2, double xi3) const {

	return (-1.0 / sqrt(2.0)) * (v1 * getYukawaMatrices_()[0] + exp(-i * xi2) * v2 * getYukawaMatrices_()[1] + exp(-i * xi3) * v3 * getYukawaMatrices_()[2]);

}

void MyVector::getMassRatio(double step, vector<double>& vd) const {

	double v1;
	Matrix3cd m;
	Matrix3d diag;
	vector<double> mume, mtme, mtmu;
	vector<vector<double>> vvd;
	vd.reserve(6);

	for (double v2 = 1.0; v2 < 246.0; v2 += step) {
		for (double v3 = v2; v3 < sqrt(pow(246.0, 2.0) - pow(v2, 2.0)); v3 += step) {
			if ((v2 == 0.0) & (v3 == 0.0)) continue;
			for (double xi2 = -1.0 * pi; xi2 < pi; xi2 += 0.5 * pi) {
				for (double xi3 = -1.0 * pi; xi3 < pi; xi3 += 0.5 * pi) {
					v1 = sqrt(pow(246.0, 2.0) - pow(v2, 2.0) - pow(v3, 2.0));
					flash(to_string(v1) + " " + to_string(v2) + " " + to_string(v3), 45, 1);
					m = getMassMatrix(v1, v2, v3, xi2, xi3); m = m * m.adjoint();
					JacobiSVD<Matrix3cd> svd(m, ComputeFullU);
					diag = svd.singularValues().asDiagonal();
					mume.emplace_back(sqrt(diag(1, 1) / diag(2, 2)));
					mtme.emplace_back(sqrt(diag(0, 0) / diag(2, 2)));
					mtmu.emplace_back(sqrt(diag(0, 0) / diag(1, 1)));
				}
			}
		}
	}
	cout << endl;

	vvd = { mume,mtme,mtmu };
	for (auto i : vvd) {
		auto minmax = minmax_element(begin(i), end(i));
		vd.insert(vd.end(), { *minmax.first,*minmax.second });
	}

}

MyVector MyVector::segment(int i, int n) const { return MyVector(v_.segment(i, n)); }

complex<double>* MyVector::array() const {

	static complex<double> array[9];
	for (int i = 0; i < 9; i++) array[i] = v_[i];
	return array;

}

vector<complex<double>> MyVector::getPhases() {

	vector<complex<double>> vcd;
	bool justonce = true;

	for (int i = 0; i < getNumberofElements(); i++) {
		if (!(abs(v_(i)) < DBL_EPSILON * 10)) {
			if (justonce) { justonce = false; continue; }
			vcd.push_back(v_(i));
		}
	}

	return vcd;

}

MyVector MyVector::setAllNonZeroElementstoPhases() const {

	vector<VectorXcd> vv;
	vector<int> vi;
	int min_index, max_index;
	VectorXcd v1(9), v2(9), v(27);

	vv = { v_.segment(0,9),v_.segment(9,9),v_.segment(18,9) };
	for (auto i : vv) {
		for (int j = 0; j < i.size(); j++) {
			if (!(abs(i(j)) < DBL_EPSILON * 10)) { vi.push_back(j); break; }
		}
	}
	min_index = min_element(vi.begin(), vi.end()) - vi.begin();
	max_index = max_element(vi.begin(), vi.end()) - vi.begin();
	v1 << vv[min_index]; v2 << vv[max_index];
	if (min_index > max_index) { vv.erase(vv.begin() + min_index); vv.erase(vv.begin() + max_index); }
	else { vv.erase(vv.begin() + max_index); vv.erase(vv.begin() + min_index); }
	v << v1, vv[0], v2; setFirstNonZeroElementto1_(v);

	return MyVector(v);

}

MyVector MyVector::setAllNonZeroElementsto1() const {

	VectorXcd v = VectorXcd::Zero(27);

	for (int i = 0; i < getNumberofElements(); i++) {
		if (!(abs(v_(i)) < DBL_EPSILON * 10)) {
			v(i).real(1.0); v(i).imag(0.0);
		}
	}

	return MyVector(v);

}

MyVector operator*(const MatrixXcd& m, const MyVector& mv) {

	return MyVector(m * mv.v_);

}

ostream& operator<<(ostream& os, const MyVector& mv) {

	return os << mv.v_.format(fmtV);

}