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

	typedef complex<double> cd;

	VectorXcd v(27);
	v << cd(0, 0), cd(0, 0), cd(0, 0), cd(0, -5.55112e-17), cd(0, 0), cd(0.0911973, -0.157958), cd(0, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(-0.182395, -5.71479e-09), cd(-5.55112e-17, -1.11022e-16), cd(0, 0),
		cd(0, 0), cd(0.0911973, 0.157958), cd(-4.16334e-17, 2.22045e-16), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0);

	MyVector mv(v);
	vector<double> vd;

	cout << "Input myvector: " << endl << mv << endl << endl;

	cout << "getNumberofElements: " << mv.getNumberofElements() << endl << endl;

	mv.setNumericZerotoActualZero();
	cout << "setNumericZerotoActualZero: " << endl << mv << endl << endl;

	mv.setFirstNonZeroElementto1();
	cout << "setFirstNonZeroElementto1: " << endl << mv << endl << endl;

	mv = mv.setAllNonZeroElementstoPhases();
	cout << "setAllNonZeroElementstoPhases: " << endl << mv << endl << endl;

	cout << "getPhases: ";  for (auto& i : mv.getPhases()) cout << i << " "; cout << endl << endl;

	cout << "getMassMatrix (with default parameters): " << endl << mv.getMassMatrix() << endl << endl;

	cout << "getMassRatio (with step equal to 10.0): " << endl;
	mv.getMassRatio(10, vd); for (auto i : vd) cout << i << " "; cout << endl << endl;

	mv = mv.setAllNonZeroElementsto1();
	cout << "setAllNonZeroElementsto1: " << endl << mv << endl;

	return EXIT_SUCCESS;

}