#pragma once
#ifndef MYMATRIX_H_
#define MYMATRIX_H_

#include "Eigen/Dense"


class MyVector;
class Yukawa;

class MyMatrix
{
public:
	MyMatrix();
	explicit MyMatrix(const Eigen::MatrixXcd&);
	int getNumberofRows() const { return row_; }
	int getNumberofCols() const { return col_; }
	void setNumericZerotoActualZero();
	bool isEigenvector1(const MyVector&) const;
	friend MyMatrix getKroneckerProduct3(const MyMatrix&, const MyMatrix&, const MyMatrix&, const std::string&);
	MyMatrix getEigenvectors1() const;
	friend MyMatrix getIntersectionBasis(const MyMatrix&, const MyMatrix&);
	Yukawa extractYukawaSolution() const;
	friend std::ostream& operator<<(std::ostream&, const MyMatrix&);
private:
	friend void setNumericZerotoActualZero_(Eigen::MatrixXcd&);
	int row_;
	int col_;
	Eigen::MatrixXcd m_;
};
MyMatrix getKroneckerProduct3(const MyMatrix&, const MyMatrix&, const MyMatrix&, const std::string&);
MyMatrix getIntersectionBasis(const MyMatrix&, const MyMatrix&);

#endif // !MYMATRIX_H_
