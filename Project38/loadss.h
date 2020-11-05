#pragma once
#ifndef loadss_h
#define loadss_h


template<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
inline bool loadMatrix(std::stringstream& ss, Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>& m, int n_rows) {

	if (ss.str() == "") {
		std::cerr << "ERROR. Cannot find stringstream!";
		m = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>(0, 0);
		return false;
	}

	Scalar d;
	std::vector<Scalar> v;
	while (ss >> d) { v.push_back(d); }

	if (v.size() % n_rows != 0) {
		std::cerr << "ERROR. Wrong number of elements!";
		m = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>(0, 0);
		return false;
	}

	int n_cols;
	n_cols = v.size() / n_rows;
	m = Eigen::Matrix<Scalar, RowsAtCompileTime, ColsAtCompileTime>(n_rows, n_cols);

	for (int i = 0; i < n_rows; i++)
		for (int j = 0; j < n_cols; j++)
			m(i, j) = v[i * n_cols + j];

	return true;

};

#endif