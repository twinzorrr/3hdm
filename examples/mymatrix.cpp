#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"

#include "basicf.h"
#include "loadm.h"
#include "myvector.h"
#include "mymatrix.h"
#include "group.h"
#include "yukawa.h"

using namespace std;
using namespace Eigen;


int main() {

	Matrix3cd m1, m2, m3;
	m1 << 0, 1, 0, 0, 0, 1, 1, 0, 0;
	m2 << -1, 0, 0, 0, 1, 0, 0, 0, -1;
	m3 << -1, 0, 0, 0, -1, 0, 0, 0, 1;
	
	MyMatrix mm1(m1), mm2(m2), mm3(m3), kp1, kp2, kp3, es1, es2, es3, ib1, ib2;

	kp1 = getKroneckerProduct3(mm1, mm1, mm1, "charged");
	kp2 = getKroneckerProduct3(mm2, mm2, mm2, "charged");
	kp3 = getKroneckerProduct3(mm3, mm3, mm3, "charged");
	cout << "getKroneckerProduct3: " << endl << kp1 << endl << endl;
	
	es1 = kp1.getEigenvectors1();
	es2 = kp2.getEigenvectors1();
	es3 = kp3.getEigenvectors1();
	cout << "getEigenvectors1: " << endl << es1 << endl << endl;
	
	es1.setNumericZerotoActualZero();
	cout << "setNumericZerotoActualZero: " << endl << es1 << endl << endl;
	
	cout << "getNumberofRows: " << es1.getNumberofRows() << " getNumberofCols: " << es1.getNumberofCols() << endl << endl;

	ib1 = getIntersectionBasis(es1, es2);
	ib2 = getIntersectionBasis(ib1, es3);
	cout << "getIntersectionBasis: " << endl << ib2 << endl << endl;
	
	cout << "isEigenvector1: ";
	cout << kp1.isEigenvector1(ib2.extractYukawaSolution()[0]) << " ";
	cout << kp2.isEigenvector1(ib2.extractYukawaSolution()[0]) << " ";
	cout << kp3.isEigenvector1(ib2.extractYukawaSolution()[0]) << endl;

	return EXIT_SUCCESS;

}