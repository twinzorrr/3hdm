#include <iostream>
#include <fstream>
#include <vector>
#include "../Eigen/Dense"

#include "../3hdm"


int main() {
	
	using cd = std::complex<double>;

	std::fstream ifile, ofile1, ofile2, ofile3;
	std::vector<std::string> vs;
	std::vector<size_t> vn;
	std::vector<Eigen::MatrixXcd> vm;
	std::vector<MyMatrix<cd>> vmm;
	std::vector<Yukawa<cd>> vyc, vyd, vyp;
	Yukawa<cd> yc, yd, yp, ypm, ypp, ypu, yppc, yppd, yppcm, yppdm;

	ofile1 = basic::fileOpener("outputs/", "charged_all.txt", std::ios::out);
	ofile2 = basic::fileOpener("outputs/", "dirac_all.txt", std::ios::out);
	ofile3 = basic::fileOpener("outputs/", "pair_all.txt", std::ios::out);

	(basic::loadString("", "groups.txt")).swap(vs);
	(basic::loadInteger("", "nors.txt")).swap(vn);

	for (size_t s = 0; s < vs.size(); s++) {
		ifile = basic::fileOpener("gs/", vs[s], std::ios::in);
		vm = load::loadM<cd>(ifile, 3, 8);
		vmm.reserve(vm.size());
		for (const auto& i : vm) vmm.emplace_back(i);
		Group<cd> g(vmm, vs[s], vn[s]); vm.clear();
		(g.findPairSolutions(vyc, vyd, ofile1, ofile2, ofile3)).swap(vyp);
		findUniqueVectors(vyc, yc, vs[s]);
		findUniqueVectors(vyd, yd, vs[s]);
		if (vyp.size()) findUniquePairs(vyp, yp, vs[s]);
		vmm.clear();
	}

	ofile3.close();
	ofile2.close();
	ofile1.close();

	yc.printToFile("charged_unique.txt");
	yd.printToFile("dirac_unique.txt");
	yp.printToFile("pair_unique.txt");
	
	findUniqueMatrixPairs(yp, ypm);
	ypm.printToFile("pair_massmatrix.txt", {});

	yp.findSolutionsWithUnityElements(ypp, ypu);
	ypp.printToFile("pair_withphases.txt");
	ypu.printToFile("pair_withunityelements.txt");
	
	ypp.splitPairsintoVectors(yppc, yppd);
	findUniqueMatrices(yppc, yppcm);
	findUniqueMatrices(yppd, yppdm);
	yppcm.printToFile("charged_massmatrix.txt", {});
	yppdm.printToFile("dirac_massmatrix.txt", {});
	
	return EXIT_SUCCESS;
	
}