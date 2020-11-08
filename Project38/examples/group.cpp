#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"

#include "../basicf.h"
#include "../loadm.h"
#include "../myvector.h"
#include "../mymatrix.h"
#include "../group.h"
#include "../yukawa.h"


int main() {

	// In general, one can obtain given group by loading appropriate file with the name given by id of this group
	// The list of all physical groups can be found in the groups.txt file, while the number of 3-dim representations
	// for these groups can be found in the nors.txt file

	std::fstream ifile;
	std::vector<Eigen::MatrixXcd> vm;
	std::vector<MyMatrix> vmm;

	ifile = fileOpener("../gs/", "[ 12, 3 ]", std::ios::in);
	vm = loadM(ifile, 3, 8);
	// the second argument determines dim of matrix while the third one determines decimal precision of its elements
	for (auto i : vm) vmm.emplace_back(i);

	Group g(vmm, "[ 12, 3 ]", 1), kp, es; // the third argument determines the number of representations of given group
	MyMatrix ib;

	// for one combination of matrix representations of given group

	kp = g.findKroneckerProduct3(0, 0, 0, "charged");
	std::cout << "findKroneckerProduct3: " << kp.getSize() << std::endl << kp[0] << std::endl << std::endl;

	es = kp.findEigenvectors1();
	std::cout << "findEigenvectors1: " << es.getSize() << std::endl << es[0] << std::endl << std::endl;

	es.setNumericZerotoActualZero();
	std::cout << "setNumericZerotoActualZero: " << std::endl << es[0] << std::endl << std::endl;

	ib = es.findIntersectionBasis();
	std::cout << "findIntersectionBasis: " << std::endl << ib << std::endl << std::endl;

	std::cout << "isEigenvector1: " << kp.isEigenvector1(ib.extractYukawaSolution()[0]) << std::endl << std::endl;

	// for all combinations of matrix representations of given group

	std::fstream ofile_c, ofile_d, ofile_ps;
	std::vector<Yukawa> vy_c, vy_d, vy_ps;

	vy_c = g.findSolutions("charged", ofile_c);
	vy_d = g.findSolutions("dirac", ofile_d);
	std::cout << "findSolutions: " << vy_c.size() << " (charged) " << vy_d.size() << " (dirac) " << std::endl << std::endl;

	// or if one also wants to find pair solutions then instead of above 2 methods we can call the following one

	vy_c.clear(); vy_d.clear();
	vy_ps = g.findPairSolutions(vy_c, vy_d, ofile_c, ofile_d, ofile_ps);
	std::cout << "findPairSolutions: " << vy_c.size() << " (charged) " << vy_d.size() << " (dirac) " << vy_ps.size() << " (pair) " << std::endl;

	// It is not recommended to initialize groups by our hands, however if one wants to do that then remember that
	// their matrix representations must be 3-dim

	return EXIT_SUCCESS;

}