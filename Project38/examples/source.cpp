#include <iostream>
#include <fstream>
#include <vector>
#include "../Eigen/Dense"

#include "../basicf.h"
#include "../loadm.h"
#include "../myvector.h"
#include "../mymatrix.h"
#include "../group.h"
#include "../yukawa.h"


int main() {
	
	std::fstream ifile, ofile1, ofile2, ofile3;
	std::vector<std::string> vs;
	std::vector<size_t> vn;
	std::vector<Eigen::MatrixXcd> vm;
	std::vector<MyMatrix> vmm;
	std::vector<Yukawa> vys_c, vys_d, vys_ps;
	Yukawa ysu_c, ysu_d, ysu_ps, ysu_psm, ys_wp, ys_wue, ys_wpc, ys_wpd, ys_wpcm, ys_wpdm;

	ofile1 = fileOpener("outputs/", "charged_all.txt", std::ios::out);
	ofile2 = fileOpener("outputs/", "dirac_all.txt", std::ios::out);
	ofile3 = fileOpener("outputs/", "pair_solutions.txt", std::ios::out);

	vs = loadString("", "groups.txt");
	vn = loadInteger("", "nors.txt");

	for (size_t s = 0; s < vs.size(); s++) {
		ifile = fileOpener("gs/", vs[s], std::ios::in);
		vm = loadM(ifile, 3, 8);
		vmm.reserve(vm.size());
		for (auto i : vm) vmm.emplace_back(i);
		Group g(vmm, vs[s], vn[s]); vm.clear();
		(g.findPairSolutions(vys_c, vys_d, ofile1, ofile2, ofile3)).swap(vys_ps);
		findUniqueVectors(vys_c, ysu_c, vs[s]);
		findUniqueVectors(vys_d, ysu_d, vs[s]);
		if (vys_ps.size()) findUniquePairs(vys_ps, ysu_ps, vs[s]);
		vmm.clear();
	}

	ofile3.close();
	ofile2.close();
	ofile1.close();

	ysu_c.printToFile("charged_unique.txt");
	ysu_d.printToFile("dirac_unique.txt");
	ysu_ps.printToFile("pair_solutions_unique.txt");
	
	findUniqueMatrixPairs(ysu_ps, ysu_psm);
	ysu_psm.printToFile("pair_mass_matrices_unique.txt", {});

	ysu_ps.findSolutionsWithUnityElements(ys_wp, ys_wue);
	ys_wp.printToFile("pair_withphases.txt");
	ys_wue.printToFile("pair_withunityelements.txt");
	
	ys_wp.splitPairsintoVectors(ys_wpc, ys_wpd);
	findUniqueMatrices(ys_wpc, ys_wpcm);
	findUniqueMatrices(ys_wpd, ys_wpdm);
	ys_wpcm.printToFile("charged_mass_matrices_unique.txt", {});
	ys_wpdm.printToFile("dirac_mass_matrices_unique.txt", {});
	
	return EXIT_SUCCESS;
	
}