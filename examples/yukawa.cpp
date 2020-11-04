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

	// In general, Yukawa vectors are obtained as a solutions to the invariant equation determined by the model we are using
	// Here findSolutions and findPairSolutions methods of group class are used to obtain vector container of yukawa objects

	fstream ifile, ofile1, ofile2, ofile3;
	vector<MatrixXcd> vm;
	vector<MyMatrix> vmm;

	ifile = fileOpener("gs/", "[ 21, 1 ]", ios::in);
	vm = loadM(ifile, 3, 8);
	for (auto i : vm) vmm.emplace_back(i);
	Group g(vmm, "[ 21, 1 ]", 2);

	vector<Yukawa> vy_c, vy_d, vy_ps;
	vy_ps = g.findPairSolutions(vy_c, vy_d, ofile1, ofile2, ofile3);
	cout << endl;

	Yukawa yu_c, yu_d;
	findUniqueVectors(vy_c, yu_c, "[ 21, 1 ]");
	findUniqueVectors(vy_c, yu_d, "[ 21, 1 ]");
	cout << "findUniqueVectors (for charged particles): " << yu_c.getSize() << endl << yu_c << endl;
	//yu_c.printToFile(filename, po::VECTOR);

	Yukawa yu_ps;
	findUniquePairs(vy_ps, yu_ps, "[ 21, 1 ]");
	cout << "findUniquePairs: " << yu_ps.getSize() / 2 << endl << yu_ps << endl;
	//yu_ps.printToFile(filename, po::PAIR);

	Yukawa yu_wp;
	yu_ps.findSolutionsWithPhases(yu_wp);
	cout << "findSolutionsWithPhases: " << yu_wp.getSize() / 2 << endl << yu_wp << endl;

	// or if one also wants to find solutions with unity element then instead of above method we can call the following one

	yu_wp.clear();
	Yukawa yu_wue;
	yu_ps.findSolutionsWithUnityElements(yu_wp, yu_wue);
	cout << "findSolutionsWithUnityElements: " << yu_wue.getSize() / 2 << endl << yu_wue << endl;

	Yukawa yu_psc, yu_psd;
	yu_ps.splitPairsintoVectors(yu_psc, yu_psd);
	cout << "splitPairsintoVectors (charged part): " << endl << yu_psc << endl;

	Yukawa yu_pscm, yu_psdm;
	findUniqueMatrices(yu_psc, yu_pscm);
	findUniqueMatrices(yu_psd, yu_psdm);
	cout << "findUniqueMatrices (charged part): " << endl << yu_pscm << endl;
	//yu_pscm.printToFile(filename);

	// or if one also wants to find mass ratio for unique matrices then instead of above 2 methods we can call the following ones

	yu_pscm.clear(); yu_psdm.clear();
	vector<vector<double>> vd_c, vd_d;
	findMassRatio(yu_psc, 20.0, yu_pscm, vd_c);
	findMassRatio(yu_psd, 20.0, yu_psdm, vd_d);
	cout << "findUniqueMatrices (charged part): " << endl << yu_pscm << endl;
	//yu_pscm.printToFile(filename, vd_c);

	// Anyway there is no reason to not initialize Yukawa vectors by our hands remembering that these vectors must be 27-dim
	// For this purpose constructor with vector<string> as a second argument is highly recommended

	typedef complex<double> cd;

	VectorXcd v1(27), v2(27), v3(27), v4(27), v5(27);
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

	vector<MyVector> vmv;
	vector<string> vs;
	vmv.emplace_back(v1); vmv.emplace_back(v2); vmv.emplace_back(v3);
	vs.emplace_back("No group"); vs.emplace_back("No group"); vs.emplace_back("No group");
	Yukawa y(vmv, vs);

	cout << "isUniqueVector: " << y.isUniqueVector(MyVector(v3), "No group") << endl;
	cout << "isUniqueVector: " << y.isUniqueVector(MyVector(v4), "No group") << endl;

	cout << "isUniquePair: " << y.isUniquePair(MyVector(v1), MyVector(v2), "No group") << endl;
	cout << "isUniquePair: " << y.isUniquePair(MyVector(v1), MyVector(v4), "No group") << endl;

	cout << "isUniqueMatrix: " << y.isUniqueMatrix(MyVector(v4), "No group") << endl;
	cout << "isUniquematrix: " << y.isUniqueMatrix(MyVector(v5), "No group") << endl;

	return EXIT_SUCCESS;

}