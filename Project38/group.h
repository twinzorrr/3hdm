#pragma once
#ifndef GROUP_H_
#define GROUP_H_

#include <vector>
#include "Eigen/Dense"


class MyVector;
class MyMatrix;
class Yukawa;

class Group
{
public:
	Group();
	explicit Group(const std::vector<MyMatrix>&, const std::string& group = "No Group", size_t nor = 0);
	size_t getSize() const { return vm_.size(); }
	void reserve(size_t i) { vm_.reserve(i); }
	void emplace_back(const Eigen::MatrixXcd& m) { vm_.emplace_back(m); }
	void push_back(const MyMatrix& mm) { vm_.push_back(mm); }
	void clear() { vm_.clear(); }
	void swap(Group& g) { (g.vm_).swap(vm_); (g.group_).swap(group_); size_t tmp; tmp = g.nor_; g.nor_ = nor_; nor_ = tmp; }
	void setNumericZerotoActualZero();
	bool isEigenvector1(const MyVector&) const;
	Group findKroneckerProduct3(size_t, size_t, size_t, const std::string&) const;
	Group findEigenvectors1() const;
	MyMatrix findIntersectionBasis() const;
	std::vector<Yukawa> findSolutions(const std::string&, std::fstream&, std::vector<std::string>* p_vs = nullptr) const;
	std::vector<Yukawa> findPairSolutions(std::vector<Yukawa>&, std::vector<Yukawa>&, std::fstream&, std::fstream&, std::fstream&) const;
	auto begin() { return vm_.begin(); }
	auto end() { return vm_.end(); }
	auto begin() const { return vm_.begin(); }
	auto end() const { return vm_.end(); }
	MyMatrix& operator[](size_t);
	const MyMatrix& operator[](size_t) const;
	friend std::ostream& operator<<(std::ostream&, const Group&);
private:
	size_t nor_;
	size_t nog_;
	std::string group_;
	std::vector<MyMatrix> vm_;
};

#endif // !GROUP_H_
