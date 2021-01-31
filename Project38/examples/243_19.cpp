#include <iostream>
#include <fstream>
#include <vector>
#include "../Eigen/Dense"

#include "../3hdm"


int main() {

	using cd = std::complex<double>;
	using cf = std::complex<float>;

	std::fstream ifile, ofile1, ofile2, ofile3;
	std::vector<Eigen::MatrixXcd> vm;
	std::vector<MyMatrix<cd>> vmm;

	ifile = basic::fileOpener("gs/", "[ 243, 19 ]", std::ios::in);
	ofile1 = basic::fileOpener("outputs/", "243_19-charged_all.txt", std::ios::out);
	ofile2 = basic::fileOpener("outputs/", "243_19-dirac_all.txt", std::ios::out);
	ofile3 = basic::fileOpener("outputs/", "243_19-pair_all.txt", std::ios::out);

	vm = load::loadM<cd>(ifile, 3, 17);
	vmm.reserve(vm.size());
	for (const auto& i : vm) vmm.emplace_back(i);
	Group<cd> g(vmm, "[ 243, 19 ]", 24);
	vmm.clear(); vm.clear();

	std::vector<Yukawa<cf>> vyc, vyd, vyp;
	Yukawa<cf> yc, yd, yp, ypp, ypu, yppc, yppd, yppm, yppcm;
	std::vector<std::vector<float>> vvf;

	(g.findPairSolutionsWithCast<cf>(vyc, vyd, ofile1, ofile2, ofile3)).swap(vyp);
	findUniqueVectors(vyc, yc, "[ 243, 19 ]");
	findUniqueVectors(vyd, yd, "[ 243, 19 ]");
	if (vyp.size()) findUniquePairs(vyp, yp, "[ 243, 19 ]");
	ofile3.close(); ofile2.close(); ofile1.close();

	yc.printToFile("243_19-charged_unique.txt");
	yd.printToFile("243_19-dirac_unique.txt");
	yp.printToFile("243_19-pair_unique.txt");

	yp.findSolutionsWithUnityElements(ypp, ypu);
	ypp.printToFile("243_19-pair_withphases.txt");
	ypu.printToFile("243_19-pair_withunityelements.txt");

	findUniqueMatrixPairs(ypp, yppm);
	yppm.printToFile("243_19-pair_massmatrix.txt", {});
	ypp.splitPairsintoVectors(yppc, yppd);
	findMassRatio(yppc, 1.0f, yppcm, vvf);
	yppcm.printToFile("243_19-charged_massratio.txt", vvf);

	return EXIT_SUCCESS;

}