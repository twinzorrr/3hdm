#pragma once
#ifndef GROUP_H_
#define GROUP_H_

#include <vector>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class MyVector;
class MyMatrix;
class Yukawa;

class Group
{
public:
	Group();
	explicit Group(vector<MyMatrix>, string group = "No Group", size_t nor = 0);
	size_t getSize() const { return vm_.size(); }
	void reserve(size_t i) { vm_.reserve(i); }
	void emplace_back(const MatrixXcd& m) { vm_.emplace_back(m); }
	void push_back(const MyMatrix& mm) { vm_.push_back(mm); }
	void swap(Group& g) { (g.vm_).swap(vm_); (g.group_).swap(group_); size_t tmp; tmp = g.nor_; g.nor_ = nor_; nor_ = tmp; }
	void setNumericZerotoActualZero();
	bool isEigenvector1(const MyVector&) const;
	Group findKroneckerProduct3(size_t, size_t, size_t, const string&) const;
	Group findEigenvectors1() const;
	MyMatrix findIntersectionBasis() const;
	vector<Yukawa> findSolutions(const string&, vector<string>&, fstream&) const;
	vector<Yukawa> findPairSolutions(vector<Yukawa>&, vector<Yukawa>&, fstream&, fstream&, fstream&) const;
	vector<MyMatrix>::iterator begin() { return vm_.begin(); }
	vector<MyMatrix>::iterator end() { return vm_.end(); }
	vector<MyMatrix>::const_iterator begin() const { return vm_.begin(); }
	vector<MyMatrix>::const_iterator end() const { return vm_.end(); }
	MyMatrix& operator[](size_t);
	const MyMatrix& operator[](size_t) const;
	friend ostream& operator<<(ostream&, const Group&);
private:
	size_t nor_;
	size_t nog_;
	string group_;
	vector<MyMatrix> vm_;
};

#endif // !GROUP_H_
