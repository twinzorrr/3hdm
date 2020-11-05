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


int main() {

	Eigen::Matrix3cd m1, m2, m3;
	m1 << 0, 1, 0, 0, 0, 1, 1, 0, 0;
	m2 << -1, 0, 0, 0, 1, 0, 0, 0, -1;
	m3 << -1, 0, 0, 0, -1, 0, 0, 0, 1;
	
	MyMatrix mm1(m1), mm2(m2), mm3(m3), kp1, kp2, kp3, es1, es2, es3, ib1, ib2;

	kp1 = getKroneckerProduct3(mm1, mm1, mm1, "charged");
	kp2 = getKroneckerProduct3(mm2, mm2, mm2, "charged");
	kp3 = getKroneckerProduct3(mm3, mm3, mm3, "charged");
	std::cout << "getKroneckerProduct3: " << std::endl << kp1 << std::endl << std::endl;
	
	es1 = kp1.getEigenvectors1();
	es2 = kp2.getEigenvectors1();
	es3 = kp3.getEigenvectors1();
	std::cout << "getEigenvectors1: " << std::endl << es1 << std::endl << std::endl;
	
	es1.setNumericZerotoActualZero();
	std::cout << "setNumericZerotoActualZero: " << std::endl << es1 << std::endl << std::endl;
	
	std::cout << "getNumberofRows: " << es1.getNumberofRows() << " getNumberofCols: " << es1.getNumberofCols() << std::endl << std::endl;

	ib1 = getIntersectionBasis(es1, es2);
	ib2 = getIntersectionBasis(ib1, es3);
	std::cout << "getIntersectionBasis: " << std::endl << ib2 << std::endl << std::endl;
	
	std::cout << "isEigenvector1: ";
	std::cout << kp1.isEigenvector1(ib2.extractYukawaSolution()[0]) << " ";
	std::cout << kp2.isEigenvector1(ib2.extractYukawaSolution()[0]) << " ";
	std::cout << kp3.isEigenvector1(ib2.extractYukawaSolution()[0]) << std::endl;

	return EXIT_SUCCESS;

}
