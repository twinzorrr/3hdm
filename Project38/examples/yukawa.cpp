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

	// In general, Yukawa vectors are obtained as a solutions to the invariant equation determined by the model we are using
	// findSolutions and findPairSolutions methods of group class are used to obtain vector container of yukawa objects
	// using findSolutions method we get Yukawa pre-vectors, while using findPairSolutions method we also get Yukawa pre-pairs

	std::fstream ifile, ofile1, ofile2, ofile3;
	std::vector<Eigen::MatrixXcd> vm;
	std::vector<MyMatrix> vmm;

	ifile = fileOpener("../gs/", "[ 21, 1 ]", std::ios::in);
	vm = loadM(ifile, 3, 8);
	// the second argument determines dim of matrix while the third one determines decimal precision of its elements
	for (auto i : vm) vmm.emplace_back(i);
	Group g(vmm, "[ 21, 1 ]", 2); // the third argument determines the number of representations of given group

	std::vector<Yukawa> vy_c, vy_d, vy_ps;
	vy_ps = g.findPairSolutions(vy_c, vy_d, ofile1, ofile2, ofile3);
	std::cout << std::endl;

	Yukawa yu_c, yu_d;
	findUniqueVectors(vy_c, yu_c, "[ 21, 1 ]");
	findUniqueVectors(vy_c, yu_d, "[ 21, 1 ]");
	std::cout << "findUniqueVectors (for charged particles): " << yu_c.getSize() << std::endl << yu_c << std::endl;
	//yu_c.printToFile(filename);

	Yukawa yu_ps;
	findUniquePairs(vy_ps, yu_ps, "[ 21, 1 ]");
	std::cout << "findUniquePairs: " << yu_ps.getSize() / 2 << std::endl << yu_ps << std::endl;
	//yu_ps.printToFile(filename);

	Yukawa yu_wp;
	yu_ps.findSolutionsWithPhases(yu_wp);
	std::cout << "findSolutionsWithPhases (for pairs): " << yu_wp.getSize() / 2 << std::endl << yu_wp << std::endl;
	
	// or if one also wants to find solutions with unity elements for pairs then instead of above method we can call the following one

	yu_wp.clear();
	Yukawa yu_wue;
	yu_ps.findSolutionsWithUnityElements(yu_wp, yu_wue);
	std::cout << "findSolutionsWithUnityElements (for pairs): " << yu_wue.getSize() / 2 << std::endl << yu_wue << std::endl;

	Yukawa yu_psm;
	findUniqueMatrixPairs(yu_ps, yu_psm);
	std::cout << "findUniqueMatrixPairs: " << yu_psm.getSize() / 2 << std::endl << yu_psm << std::endl;
	// yu_psm.printToFile(filename, {});

	// or if one also wants to find mass ratio for unique pairs of matrices then instead of above method we can call the following one

	yu_psm.clear();
	std::vector<std::vector<double>> vd_p;
	findMassRatio(yu_ps, 20.0, yu_psm, vd_p);
	std::cout << "findMassRatio (for pairs): " << std::endl << yu_psm << std::endl;
	// yu_psm.printToFile(filename, vd_p);

	Yukawa yu_psc, yu_psd;
	yu_ps.splitPairsintoVectors(yu_psc, yu_psd);
	std::cout << "splitPairsintoVectors (charged vectors part): " << std::endl << yu_psc << std::endl;

	// one can also find solutions with phases and unity elements as well unique matrices for Yukawa vectors as follows

	Yukawa yu_wpc, yu_wpd;
	yu_psc.findSolutionsWithPhases(yu_wpc);
	yu_psd.findSolutionsWithPhases(yu_wpd);
	std::cout << "findSolutionsWithPhases (for charged vectors): " << yu_wpc.getSize() << std::endl << yu_wpc << std::endl;

	// or if one also wants to find solutions with unity elements for vectors then instead of above 2 methods we can call the following ones

	yu_wpc.clear(); yu_wpd.clear();
	Yukawa yu_wuec, yu_wued;
	yu_psc.findSolutionsWithUnityElements(yu_wpc, yu_wuec);
	yu_psd.findSolutionsWithUnityElements(yu_wpd, yu_wued);
	std::cout << "findSolutionsWithUnityElements (for charged vectors): " << yu_wpc.getSize() << std::endl << yu_wpc << std::endl;

	Yukawa yu_pscm, yu_psdm;
	findUniqueMatrices(yu_psc, yu_pscm);
	findUniqueMatrices(yu_psd, yu_psdm);
	std::cout << "findUniqueMatrices (charged part): " << std::endl << yu_pscm << std::endl;
	//yu_pscm.printToFile(filename, {});

	// or if one also wants to find mass ratio for unique matrices then instead of above 2 methods we can call the following ones

	yu_pscm.clear(); yu_psdm.clear();
	std::vector<std::vector<double>> vd_c, vd_d;
	findMassRatio(yu_psc, 20.0, yu_pscm, vd_c);
	findMassRatio(yu_psd, 20.0, yu_psdm, vd_d);
	std::cout << "findMassRatio (charged part): " << std::endl << yu_pscm << std::endl;
	//yu_pscm.printToFile(filename, vd_c);

	// Anyway there is no reason to not initialize Yukawa vectors by our hands remembering that these vectors must be 27-dim
	// For this purpose constructor with yo, vector<MyVector> and vector<string> as a arguments is absolutely necessary
	// yo determines type of solution and can take one of two values, yo::VECTOR (for vector solutions) or yo::PAIR (for par solutions)

	using cd = std::complex<double>;

	Eigen::VectorXcd v1(27), v2(27), v3(27), v4(27), v5(27);
	v1 << cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0);
	v2 << cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0);
	v3 << cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(-0.5, 0.866025), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0);
	v4 << cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(-0.5, -0.866025), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0);
	v5 << cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(1, 0), cd(0, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(-0.5, -0.866025), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(-0.5, 0.866025), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0);

	std::vector<MyVector> vmv;
	std::vector<std::string> vs;
	vmv.emplace_back(v1); vmv.emplace_back(v2); vmv.emplace_back(v3);
	vs.emplace_back("No group"); vs.emplace_back("No group"); vs.emplace_back("No group");
	Yukawa yv(yo::VECTOR, vmv, vs), yp(yo::PAIR, vmv, vs);

	std::cout << "isUniqueVector: " << yv.isUniqueVector(MyVector(v3), "No group") << std::endl;
	std::cout << "isUniqueVector: " << yv.isUniqueVector(MyVector(v4), "No group") << std::endl;

	std::cout << "isUniquePair: " << yp.isUniquePair(MyVector(v1), MyVector(v2), "No group") << std::endl;
	std::cout << "isUniquePair: " << yp.isUniquePair(MyVector(v1), MyVector(v4), "No group") << std::endl;

	std::cout << "isUniqueMatrix: " << yv.isUniqueMatrix(MyVector(v4), "No group") << std::endl;
	std::cout << "isUniqueMatrix: " << yv.isUniqueMatrix(MyVector(v5), "No group") << std::endl;

	std::cout << "isUniqueMatrixPair: " << yp.isUniqueMatrixPair(MyVector(v1), MyVector(v2), "No group") << std::endl;
	std::cout << "isUniqueMatrixPair: " << yp.isUniqueMatrixPair(MyVector(v4), MyVector(v4), "No group") << std::endl;

	return EXIT_SUCCESS;

}