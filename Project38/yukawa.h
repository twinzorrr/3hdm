#pragma once
#ifndef YUKAWA_H_
#define YUKAWA_H_

#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class MyVector;
class Group;

using po = enum class print_option { VECTOR = 1, PAIR };


class Yukawa
{
public:
	Yukawa();
	explicit Yukawa(vector<MyVector>, vector<string> vs = {});
	size_t getSize() const { return vv_.size(); }
	void reserve(size_t i) { vv_.reserve(i); }
	void emplace_back(const VectorXcd& v) { vv_.emplace_back(v); }
	void push_back(const MyVector& mv) { vv_.push_back(mv); }
	void clear() { vv_.clear(); vs_.clear(); }
	void swap(Yukawa& y) { y.vv_.swap(vv_); y.vs_.swap(vs_); }
	void setNumericZerotoActualZero();
	void setFirstNonZeroElementto1();
	bool errorGuard(const Group&) const;
	void splitPairsintoVectors(Yukawa&, Yukawa&) const;
	bool isUniqueVector(const MyVector&, const string& s = "");
	friend void findUniqueVectors(vector<Yukawa>&, Yukawa&, const string&);
	bool isUniquePair(const MyVector&, const MyVector&, const string& s = "");
	friend void findUniquePairs(vector<Yukawa>&, Yukawa&, const string&);
	bool isUniqueMatrix(const MyVector&, const string& s = "");
	friend void findUniqueMatrices(const Yukawa&, Yukawa&);
	bool isUniqueMatrixPair(const MyVector&, const MyVector&, const string&);
	friend void findUniqueMatrixPairs(const Yukawa&, Yukawa&);
	friend void findMassRatio(const Yukawa&, double, Yukawa&, vector<vector<double>>&, po);
	void findSolutionsWithPhases(Yukawa&, po) const;
	void findSolutionsWithUnityElements(Yukawa&, Yukawa&, po) const;
	void printToFile(const string&, po);
	void printToFile(const string&, po, vector<vector<double>> vvd);
	auto begin() { return vv_.begin(); }
	auto end() { return vv_.end(); }
	auto begin() const { return vv_.begin(); }
	auto end() const { return vv_.end(); }
	MyVector& operator[](size_t);
	const MyVector& operator[](size_t) const;
	friend ostream& operator<<(ostream&, const Yukawa&);
private:
	vector<Matrix3cd> convertSolutionstoMassMatrices() const;
	void printMassRatio(vector<double>&, fstream&) const;
	vector<string> vs_;
	vector<MyVector> vv_;
};
void findUniqueVectors(vector<Yukawa>&, Yukawa&, const string&);
void findUniquePairs(vector<Yukawa>&, Yukawa&, const string&);
void findUniqueMatrices(const Yukawa&, Yukawa&);
void findUniqueMatrixPairs(const Yukawa&, Yukawa&);
void findMassRatio(const Yukawa&, double, Yukawa&, vector<vector<double>>&, po);

#endif // !YUKAWA_H_
