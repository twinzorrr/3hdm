#pragma once
#ifndef LOADSS_H_
#define LOADSS_H_

template<typename T> using Real = typename Eigen::NumTraits<T>::Real;


namespace load {
	namespace detail {

		std::regex re() {

			std::string s1, s2, s3, s4, s5, s6, s7;
			s1 = "\\s?-?(\\d+)/(\\d+)\\*E\\((\\d+)\\)\\^(\\d+)";
			s2 = "\\s?-?(\\d+)/(\\d+)\\*E\\((\\d+)\\)";
			s3 = "\\s?-?(\\d+)\\*E\\((\\d+)\\)\\^(\\d+)";
			s4 = "\\s?-?(\\d+)\\*E\\((\\d+)\\)";
			s5 = "\\s-?(\\d+)";
			s6 = ".?E\\((\\d+)\\)\\^(\\d+)";
			s7 = ".?E\\((\\d+)\\)";
			std::regex s(s1 + "|" + s2 + "|" + s3 + "|" + s4 + "|" + s5 + "|" + s6 + "|" + s7);

			return s;

		}

		template<typename _real>
		std::vector<_real> d(const std::smatch& res) {

			using boost::lexical_cast;

			int i = 0;
			_real d[4];

			d[1] = (_real)1; d[2] = (_real)1; d[3] = (_real)1;
			for (size_t r = 1; r < res.size(); r++) {
				if (res[r] != "") { d[i] = lexical_cast<_real>(res[r]); i++; }
			}

			return { d[0],d[1],d[2],d[3] };

		}

		template<typename _scalar>
		_scalar e(const std::string& s, Real<_scalar> d1, Real<_scalar> d2, Real<_scalar> d3, Real<_scalar> d4) {

			consts<Real<_scalar>> cs;
			_scalar c, d;
			bool con;

			d = { (Real<_scalar>)2 * d2 / d1, (Real<_scalar>)0 }; c = { (Real<_scalar>)1, (Real<_scalar>)0 };
			if (s.find("/") != std::string::npos) { d = { (Real<_scalar>)2 * d4 / d3, (Real<_scalar>)0 }; c = { d1 / d2, (Real<_scalar>)0 }; }
			else if (s.find("*") != std::string::npos) { d = { (Real<_scalar>)2 * d3 / d2, (Real<_scalar>)0 }; c = { d1, (Real<_scalar>)0 }; }

			con = (s[0] == '-') || (s[1] == '-');
			if (con) return -(c * exp(cs.i * cs.pi * d));
			return c * exp(cs.i * cs.pi * d);

		}

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

	}
}

#endif // !LOADSS_H_