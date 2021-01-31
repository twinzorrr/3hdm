#pragma once
#ifndef CONSTANTS_H_
#define CONSTANTS_H_


const Eigen::IOFormat fmtV(6, Eigen::DontAlignCols, " ", " ", "", "", "[", "]");
const Eigen::IOFormat fmtM(6, 0, " ", "\n", "[", "]");

template<typename T>
struct consts
{
	const std::complex<T> i = { (T)0,(T)1 };
	const T pi = acos((T)-1);
};

#endif // !CONSTANTS_H_