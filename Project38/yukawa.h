#pragma once
#ifndef YUKAWA_H_
#define YUKAWA_H_

#include <vector>
#include "Eigen/Dense"


class MyVector;
class Group;

using yo = enum class yukawa_option { UNDEFINED, VECTOR, PAIR };

class Yukawa
{
public:
	Yukawa();
	Yukawa(yo, std::vector<MyVector>& vmv, const std::vector<std::string>& vs = {});
	size_t getSize() const { return vv_.size(); }
	void reserve(size_t i) { vv_.reserve(i); }
	void emplace_back(const Eigen::VectorXcd& v) { vv_.emplace_back(v); }
	void push_back(const MyVector& mv) { vv_.push_back(mv); }
	void clear() { vv_.clear(); vs_.clear(); }
	void swap(Yukawa& y) { y.vv_.swap(vv_); y.vs_.swap(vs_); }
	void setNumericZerotoActualZero();
	void setFirstNonZeroElementto1();
	bool errorGuard(const Group&) const;
	void splitPairsintoVectors(Yukawa&, Yukawa&) const;
	bool isUniqueVector(const MyVector&, const std::string& s = "");
	friend void findUniqueVectors(std::vector<Yukawa>&, Yukawa&, const std::string&);
	bool isUniquePair(const MyVector&, const MyVector&, const std::string& s = "");
	friend void findUniquePairs(std::vector<Yukawa>&, Yukawa&, const std::string&);
	bool isUniqueMatrix(const MyVector&, const std::string& s = "");
	friend void findUniqueMatrices(const Yukawa&, Yukawa&);
	bool isUniqueMatrixPair(const MyVector&, const MyVector&, const std::string&);
	friend void findUniqueMatrixPairs(const Yukawa&, Yukawa&);
	friend void findMassRatio(const Yukawa&, double, Yukawa&, std::vector<std::vector<double>>&);
	void findSolutionsWithPhases(Yukawa&) const;
	void findSolutionsWithUnityElements(Yukawa&, Yukawa&) const;
	void printToFile(const std::string&);
	void printToFile(const std::string&, std::vector<std::vector<double>> vvd);
	auto begin() { return vv_.begin(); }
	auto end() { return vv_.end(); }
	auto begin() const { return vv_.begin(); }
	auto end() const { return vv_.end(); }
	MyVector& operator[](size_t);
	const MyVector& operator[](size_t) const;
	friend std::ostream& operator<<(std::ostream&, const Yukawa&);
private:
	std::vector<Eigen::Matrix3cd> convertSolutionstoMassMatrices() const;
	void printMassRatio(std::vector<double>&, std::fstream&) const;
	yo solution_;
	std::vector<std::string> vs_;
	std::vector<MyVector> vv_;
};
void findUniqueVectors(std::vector<Yukawa>&, Yukawa&, const std::string&);
void findUniquePairs(std::vector<Yukawa>&, Yukawa&, const std::string&);
void findUniqueMatrices(const Yukawa&, Yukawa&);
void findUniqueMatrixPairs(const Yukawa&, Yukawa&);
void findMassRatio(const Yukawa&, double, Yukawa&, std::vector<std::vector<double>>&);

#endif // !YUKAWA_H_
