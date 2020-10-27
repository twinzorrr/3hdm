#pragma once
#ifndef YUKAWA_H_
#define YUKAWA_H_

#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class MyVector;
class Group;

typedef enum class print_option { VECTOR = 1, PAIR } po;


class Yukawa
{
public:
	Yukawa();
	explicit Yukawa(vector<MyVector>);
	size_t getSize() const { return vv_.size(); }
	void reserve(size_t i) { vv_.reserve(i); }
	void emplace_back(const VectorXcd& v) { vv_.emplace_back(v); }
	void push_back(const MyVector& mv) { vv_.push_back(mv); }
	void swap(Yukawa& y) { y.vv_.swap(vv_); y.vs_.swap(vs_); }
	void setNumericZerotoActualZero();
	void setFirstNonZeroElementto1();
	bool errorGuard(const Group&) const;
	void splitPairsintoVectors(Yukawa&, Yukawa&) const;
	void uniqueVector(const string& s = "");
	friend void findUniqueVectors(vector<Yukawa>&, Yukawa&, const string&);
	void uniquePair(const string& s = "");
	friend void findUniquePairs(vector<Yukawa>&, Yukawa&, const string&);
	void uniqueMatrix(const string& s = "");
	friend void findUniqueMatrices(const Yukawa&, Yukawa&);
	friend void findMassRatio(const Yukawa&, double, Yukawa&, vector<vector<double>>&);
	void findSolutionsWithPhases(Yukawa&) const;
	void findSolutionsWithUnityElements(Yukawa&, Yukawa&) const;
	void printToFile(const string&, po);
	void printToFile(const string&, vector<vector<double>> vvd = {});
	vector<MyVector>::iterator begin() { return vv_.begin(); }
	vector<MyVector>::iterator end() { return vv_.end(); }
	vector<MyVector>::const_iterator begin() const { return vv_.begin(); }
	vector<MyVector>::const_iterator end() const { return vv_.end(); }
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
void findMassRatio(const Yukawa&, double, Yukawa&, vector<vector<double>>&);

#endif // !YUKAWA_H_
