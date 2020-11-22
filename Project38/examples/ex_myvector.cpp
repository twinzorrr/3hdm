#include <iostream>
#include <vector>
#include "../Eigen/Dense"

#include "../3hdm"


int main() {

	using cd = std::complex<double>;
	using cf = std::complex<float>;

	Eigen::VectorXcf v_cf(27);
	v_cf << cf(0, 0), cf(0, 0), cf(0, 0), cf(0, -5.55112e-08f), cf(0, 0), cf(0.0911973f, -0.157958f), cf(0, 0), cf(0, 0), cf(0, 0),
		cf(0, 0), cf(0, 0), cf(0, 0), cf(0, 0), cf(0, 0), cf(0, 0), cf(-0.182395f, -5.71479e-08f), cf(-5.55112e-08f, -1.11022e-08f), cf(0, 0),
		cf(0, 0), cf(0.0911973f, 0.157958f), cf(-4.16334e-08f, 2.22045e-08f), cf(0, 0), cf(0, 0), cf(0, 0), cf(0, 0), cf(0, 0), cf(0, 0);

	Eigen::VectorXcd v_cd(27);
	v_cd << cd(0, 0), cd(0, 0), cd(0, 0), cd(0, -5.55112e-17), cd(0, 0), cd(0.0911973, -0.157958), cd(0, 0), cd(0, 0), cd(0, 0),
		cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(-0.182395, -5.71479e-09), cd(-5.55112e-17, -1.11022e-16), cd(0, 0),
		cd(0, 0), cd(0.0911973, 0.157958), cd(-4.16334e-17, 2.22045e-16), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0), cd(0, 0);

	MyVector<std::complex<float>> mv(v_cf);
	std::vector<float> vf;

	std::cout << "Input myvector: " << std::endl << mv << std::endl << std::endl;

	std::cout << "getNumberofElements: " << mv.getNumberofElements() << std::endl << std::endl;

	mv.setNumericZerotoActualZero();
	std::cout << "setNumericZerotoActualZero: " << std::endl << mv << std::endl << std::endl;

	mv.setFirstNonZeroElementto1();
	std::cout << "setFirstNonZeroElementto1: " << std::endl << mv << std::endl << std::endl;

	mv = mv.setAllNonZeroElementstoPhases();
	std::cout << "setAllNonZeroElementstoPhases: " << std::endl << mv << std::endl << std::endl;

	std::cout << "getPhases: ";  for (auto& i : mv.getPhases()) std::cout << i << " "; std::cout << std::endl << std::endl;

	std::cout << "getMassMatrix (with default parameters): " << std::endl << mv.getMassMatrix() << std::endl << std::endl;

	std::cout << "getMassRatio (with step equal to 10.0): " << std::endl;
	mv.getMassRatio(10.0f, vf); for (auto i : vf) std::cout << i << " "; std::cout << std::endl << std::endl;

	mv = mv.setAllNonZeroElementsto1();
	std::cout << "setAllNonZeroElementsto1: " << std::endl << mv << std::endl;
	
	return EXIT_SUCCESS;

}