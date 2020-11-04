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

	// In general, one can obtain given group by loading appropriate file with the name given by id of this group
	// The list of all physical groups can be found in the groups.txt file, while the number of 3-dim representations
	// for these groups can be found in the nors.txt file

	fstream ifile;
	vector<MatrixXcd> vm;
	vector<MyMatrix> vmm;

	ifile = fileOpener("gs/", "[ 12, 3 ]", ios::in);
	vm = loadM(ifile, 3, 8);
	for (auto i : vm) vmm.emplace_back(i);

	Group g(vmm, "[ 12, 3 ]", 1), kp, es; // the third argument determines the number of representations of given group
	MyMatrix ib;

	// for one combination of matrix representations of given group

	kp = g.findKroneckerProduct3(0, 0, 0, "charged");
	cout << "findKroneckerProduct3: " << kp.getSize() << endl << kp[0] << endl << endl;

	es = kp.findEigenvectors1();
	cout << "findEigenvectors1: " << es.getSize() << endl << es[0] << endl << endl;

	es.setNumericZerotoActualZero();
	cout << "setNumericZerotoActualZero: " << endl << es[0] << endl << endl;

	ib = es.findIntersectionBasis();
	cout << "findIntersectionBasis: " << endl << ib << endl << endl;

	cout << "isEigenvector1: " << kp.isEigenvector1(ib.extractYukawaSolution()[0]) << endl << endl;

	// for all combinations of matrix representations of given group

	fstream ofile_c, ofile_d, ofile_ps;
	vector<Yukawa> vy_c, vy_d, vy_ps;

	vy_c = g.findSolutions("charged", ofile_c);
	vy_d = g.findSolutions("dirac", ofile_d);
	cout << "findSolutions: " << vy_c.size() << " (charged) " << vy_d.size() << " (dirac) " << endl << endl;

	// or if one also wants to find pair solutions then instead of above 2 methods we can call the following one
	
	vy_c.clear(); vy_d.clear();
	vy_ps = g.findPairSolutions(vy_c, vy_d, ofile_c, ofile_d, ofile_ps);
	cout << "findPairSolutions: " << vy_c.size() << " (charged) " << vy_d.size() << " (dirac) " << vy_ps.size() << " (pair) " << endl;

	// It is not recommended to initialize groups by our hands, however if one wants to do that then remember that
	// their matrix representations must be 3-dim

	return EXIT_SUCCESS;

}