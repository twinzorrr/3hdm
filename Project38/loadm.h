#pragma once
#ifndef loadm_h
#define loadm_h

#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;


vector<MatrixXcd> loadM(fstream&, int, int);

#endif