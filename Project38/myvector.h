#pragma once
#ifndef MYVECTOR_H_
#define MYVECTOR_H_

#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;


class MyVector
{
public:
	MyVector();
	explicit MyVector(VectorXcd);
	int getNumberofElements() const { return size_; }
	void setNumericZerotoActualZero();
	void setFirstNonZeroElementto1();
	double norm() const { return v_.norm(); }
	bool isApprox(const MyVector&, double) const;
	Matrix3cd getMassMatrix(double v1 = 1.0, double v2 = 1.0, double v3 = 1.0, double xi2 = 0.5 * acos(-1.0), double xi3 = 0.5 * acos(-1.0)) const;
	void getMassRatio(double, vector<double>&) const;
	MyVector segment(int, int) const;
	complex<double>* array() const;
	vector<complex<double>> getPhases();
	MyVector setAllNonZeroElementstoPhases() const;
	MyVector setAllNonZeroElementsto1() const;
	friend MyVector operator*(const MatrixXcd&, const MyVector&);
	friend ostream& operator<<(ostream&, const MyVector&);
private:
	friend void setNumericZerotoActualZero_(VectorXcd&);
	friend void setFirstNonZeroElementto1_(VectorXcd&);
	vector<Matrix3cd> getYukawaMatrices_() const;
	int size_;
	VectorXcd v_;
};

#endif // !MYVECTOR_H_
