#pragma once
#ifndef MYVECTOR_H_
#define MYVECTOR_H_

#include <vector>
#include "Eigen/Dense"


class MyVector
{
public:
	MyVector();
	explicit MyVector(const Eigen::VectorXcd&);
	int getNumberofElements() const { return size_; }
	void setNumericZerotoActualZero();
	void setFirstNonZeroElementto1();
	double norm() const { return v_.norm(); }
	bool isApprox(const MyVector&, double) const;
	Eigen::Matrix3cd getMassMatrix(double v1 = 1.0, double v2 = 1.0, double v3 = 1.0, double xi2 = 0.5 * acos(-1.0), double xi3 = 0.5 * acos(-1.0)) const;
	void getMassRatio(double, std::vector<double>&) const;
	std::vector<std::complex<double>> getPhases();
	MyVector setAllNonZeroElementstoPhases() const;
	MyVector setAllNonZeroElementsto1() const;
	friend MyVector operator*(const Eigen::MatrixXcd&, const MyVector&);
	friend std::ostream& operator<<(std::ostream&, const MyVector&);
private:
	friend int getIndexofFirstNonZeroElement(const Eigen::VectorXcd&);
	friend void setNumericZerotoActualZero_(Eigen::VectorXcd&);
	friend void setFirstNonZeroElementto1_(Eigen::VectorXcd&);
	std::vector<Eigen::Matrix3cd> getYukawaMatrices_() const;
	int size_;
	Eigen::VectorXcd v_;
};

#endif // !MYVECTOR_H_
