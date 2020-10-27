#pragma once
#ifndef MYMATRIX_H_
#define MYMATRIX_H_

#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class MyVector;
class Yukawa;


class MyMatrix
{
public:
	MyMatrix();
	explicit MyMatrix(MatrixXcd);
	int getNumberofRows() const { return row_; }
	int getNumberofCols() const { return col_; }
	void setNumericZerotoActualZero();
	bool isEigenvector1(const MyVector&) const;
	friend MyMatrix getKroneckerProduct3(const MyMatrix&, const MyMatrix&, const MyMatrix&, const string&);
	MyMatrix getEigenvectors1() const;
	friend MyMatrix getIntersectionBasis(const MyMatrix&, const MyMatrix&);
	Yukawa extractYukawaSolution() const;
	friend ostream& operator<<(ostream&, const MyMatrix&);
private:
	friend void setNumericZerotoActualZero_(MatrixXcd&);
	int row_;
	int col_;
	MatrixXcd m_;
};

#endif // !MYMATRIX_H_
