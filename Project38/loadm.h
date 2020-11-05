#pragma once
#ifndef loadm_h
#define loadm_h

#include <vector>
#include "Eigen/Dense"


std::vector<Eigen::MatrixXcd> loadM(std::fstream&, int, int);

#endif