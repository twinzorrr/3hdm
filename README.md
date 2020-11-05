# 3hdm
Yukawa solutions seeker for 3 Higgs doublet model based on idea of flavor symmetry.

To get more information about the physics behind this project read the document "3hdm.pdf".

The simplest and the best way to run the program is the "Open with visual studio" option under the green Code button in the upper right corner.

A few words about the program itself: four classes were introduced
- myvector - (at this moment) represents VectorXcd from eigen library<br/>
 > **getNumberofElements** - returns the number of elements of a given myvector<br/>
 > **setNumericZerotoActualZero** - sets numeric zeros to actual zeros, mainly used before printing myvector to file<br/>
 > **setFirstNonZeroElementto1** - multiplies myvector by scalar in such a way that the first non-zero element is equal to 1<br/>
 > **norm** - returns a norm of a given myvector<br/>
 > **isApprox** - compares a given myvector with another passed as the first argument for equality with approximation given by the second argument<br/>
 > **getMassMatrix** - for a given myvector creates a mass matrix whose form depends on the physics model, model parameters constitute the arguments list<br/>
 > **getMassRatio** - for a given myvector calculates mass ratio and returns the result in a second argument, the first argument determines the step used in calculation<br/>
 > **segment** - returns myvector containing the elements of initial myvector whose number is determined by the first argument, starting from index determined by the second argument<br/>
 > **array** - converts myvector into an array and returns a pointer to this array<br/>
 > **getPhases** - for a given myvector returns all non-zero elements except the first one<br/>
 > **setAllNonZeroElementstoPhases** - arranges submyvectors of a given myvector in ascending order with respect to the index of their first non-zero element and then multiplies myvector by a scalar to set first non-zero element to 1<br/>
 > **setAllNonZeroElementsto1** - sets all non-zero elements to 1

- mymatrix - (at this moment) represents MatrixXcd from eigen library
 > **getNumberofRows** - returns the number of rows of a given mymatrix<br/>
 > **getNumberofCols** - returns the number of columns of a given mymatrix<br/>
 > **setNumericZerotoActualZero** - sets numeric zeros to actual zeros, mainly used before printing mymatrix to file<br/>
 > **isEigenvector1** - checks if myvector passed as an argument is the eigenvector corresponding to the eigenvalue 1 of a given mymatrix<br/>
 > **getKroneckerProduct3** - returns mymatrix as a result of kronecker product of three mymatrices passed as a first three arguments, the 4th argument determines the final form of this kronecker product and may take one of two values, "charged" or "dirac"<br/>
 > **getEigenvectors1** - returns mymatrix whose columns are eigenvectors corresponding to the eigenvalue 1 of a given mymatrix<br/>
 > **getIntersectionBasis** - returns mymatrix whose columns are basis vectors of the common eigenspace of two subspaces represented by mymatrices passed as arguments<br/>
 > **extractYukawaSolution** - extracts Yukawa solutions from mymatrix obtained, in fact, as a result of finding basis for the intersection of two or more subspaces

- group - represents a vector container of mymatrices, from physical point of view represents the discrete Group itself (each element of vector container is a matrix representation of a given Group generator, but not always)
 > **getSize** - returns the number equal to the product of number of generators and matrix representations of a given Group<br/>
 > **reserve** - reserves the memory for the number of group elements given as a argument<br/>
 > **emplace_back** - adds a new element given as a argument at the end of the group without triggering copy constructor<br/>
 > **push_back** - adds a new element given as a argument at the end of the group<br/>
 > **clear** - remove all the elements of a given group, number of reps and group id remain untouched<br/>
 > **swap** - swaps the elements of two groups<br/>
 > **setNumericZerotoActualZero** - sets numeric zeros to actual zeros for all group mymatrix elements, mainly used before printing group to file<br/>
 > **isEigenvector1** - checks if myvector passed as an argument is the eigenvector corresponding to the eigenvalue 1 for each matrix of given combination of matrix representations<br/>
 > **findKroneckerProduct3** - returns group as a result of kronecker product of three mymatrices for a given combination of matrix representations determined by first three arguments, the fourth argument determines the final form of this kronecker product and may take one of two values, "charged" or "dirac"<br/>
 > **findEigenvectors1** - returns group with matrices whose columns are eigenvectors corresponding to the eigenvalue 1 for each mymatrix of a given group, in fact calculated for a given combination of matrix representations<br/>
 > **findIntersectionBasis** - returns mymatrix whose columns are basis vectors of the common eigenspace of each individual subspace represented by mymatrices of a given group, in fact calculated for a given combination of matrix representations<br/>
 > **findSolutions** - for a given Group returns Yukawa solutions for all possible combinations of matrix representations, for charged lepton and dirac neutrino particles separately, the first argument can take one of two values, "charged" or "dirac", and determines for which particles the calculations are performed, solutions are stored in the vector of yukawa objects and printed to the file passed as a second argument<br/>
 > **findPairSolutions** - for a given group returns pairs of lepton charged - dirac neutrino Yukawa solutions for all possible combinations of matrix representations, solutions are stored in the vector of yukawa objects and printed to the file passed as a fifth argument, also passes as the first and second argument a Yukawa solutions for charged lepton and dirac nuetrino particles as well printed these solutions to the file passed as a third and fourth argument, respectively

- yukawa - represents a vector container of myvectors, from physical point of view represents the vector of Yukawa solutions itself
 > **getSize** - returns the number of Yukawa solutions obtained for a given group and unique Yukawa solutions for all considered groups<br/>
 > **reserve** - reserves the memory for the number of Yukawa solutions given as a argument<br/>
 > **emplace_back** - adds a new element given as a argument at the end of the yukawa without triggering copy constructor<br/>
 > **push_back** - adds a new element given as argument at the end of the yukawa<br/>
 > **clear** - remove all the elements of a given yukawa and the group id information for these elements<br/>
 > **swap** - swaps the elements of two yukawas<br/>
 > **setNumericZerotoActualZero** - sets numeric zeros to actual zeros for all Yukawa solutions, mainly used before printing yukawa to file<br/>
 > **setFirstNonZeroElementto1** - multiplies each Yukawa solutions by scalar in such a way that the first non-zero element is equal to 1<br/>
 > **errorGuard** - checks if obtained vector solution is actual Yukawa solution in three ways: is a eigenvector for the eigenvalue 1 for a kronecker product represented by group and passed as a argument, is a zero vector, is a duplicated solution; if it's not returns false, otherwise returns true<br/>
 > **splitPairsintoVectors** - split pairs of Yukawa solutions into two Yukawa solutions passed as a arguments<br/>
 > **isUniqueVector** - checks if myvector given as a first argument is the unique one with compare to the elements of a given yukawa, if it's unique returns true and adds this myvector to yukawa, otherwise return false and does nothing<br/>
 > **findUniqueVectors** - filters the contents of vector of yukawa objects passed as the first argument and obtained for a given group and return unique Yukawa solutions stored in a yukawa object in second argument, the id of the group is given as a third argument<br/>
 > **isUniquePair** - checks if myvectors given as a first and second arguments form a unique pair with compare to all the two consecutive elements of given yukawa, if it's unique returns true and adds these myvectors to yukawa, otherwise return false and does nothing<br/>
 > **findUniquePairs** - filters the contents of vector of yukawa objects passed as a first argument and obtained for a given group and returns unique pairs of Yukawa solutions stored in a yukawa object in second argument, the id of the group is given as a third argument<br/>
 > **isUniqueMatrix** - checks if myvector given as a first argument forms the unique squared mass matrix with compare to the elements of a given yukawa, if corresponding matrix is unique returns true and adds this myvector to yukawa, otherwise return false and does nothing<br/>
 > **findUniqueMatrices** - filters the contents of yukawa passed as a first argument and returns Yukawa solutions stored in yukawa object in the second argument that form unique squared mass matrices<br/>
 > **isUniqueMatrixPair** - checks if myvectors given as a first and second arguments form a a unique pair of squared mass matrices with compare to all the two consecutive elements of given yukawa, if it's unique returns true and adds these myvectors to yukawa, otherwise return false and does nothing<br/>
 > **findUniqueMatrixPairs** - filters the contents of yukawa passed as a first argument and returns Yukawa solutions stored in yukawa object in the second argument that form unique pairs of squared mass matrices<br/>
 > **findMassRatio** - filters the contents of yukawa passed as a first argument, returns Yukawa solutions stored in yukawa object in the third argument that form unique squared mass matrices or pairs of unique squared mass matrices and calculates for them mass ratio that are returned in the fourth argument, the second argument determines the step used in the mass ratio calculation<br/>
 > **findSolutionsWithPhases** - for a given yukawa searches for unique Yukawa solutions or unique pairs of Yukawa solutions prepared in such a way that the appropriate submyvectors are arranged in ascending order with respect to the index of their first non-zero element and with the first non-zero element of such modified myvector equal to 1, the results are returned in the first argument<br/>
 > **findSolutionsWithUnityElements** - for a given yukawa searches for unique Yukawa solutions or unique pairs of Yukawa solutions prepared in such a way that the appropriate submyvectors are arranged in ascending order with respect to the index of their first non-zero element and with the all non-zero elements of such modified myvector equal to 1, the results are returned in the second argument while as a first argument unique solutions with only first non-zero element equal to 1 are returned<br/>
 > **printToFile** - for a given yukawa prints Yukawa solutions or pairs of Yukawa solutions to the file which name is defined in the first argument<br/>
 > **printToFile** - for a given yukawa prints squared mass matrices or pairs of squared mass matrices to the file which name is defined in the first argument with or without mass ratio depends on the second argument, if one call this method and initialize the second argument with empty list, {}, then only squared mass matrices will be printed
