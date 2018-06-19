#ifndef MATRIX_OPERATIONS_H
#define MATRIX_OPERATIONS_H
#include <vector>
#include <algorithm>
#include <iostream>
class MXOP
{
public:
	static std::vector<double> triu(const std::vector<std::vector<double> >& A);
	static std::vector<double> diag(const std::vector<std::vector<double> >& A);
	static bool size(unsigned int&n, unsigned int&m, const std::vector<std::vector<double> >& matrix);
	static std::vector<std::vector<double> > zeros(const unsigned int& n, const unsigned int& m);
	static void printMatrix(const std::vector<std::vector<double> >& matrix);
	static void printVector(const std::vector<double>&vector);
	static std::vector<std::vector<double> > multiply(const std::vector<std::vector<double> >& A, const std::vector<std::vector<double> >& B);
	static std::vector<std::vector<double> > transpose(const std::vector<std::vector<double> >& mat);
	static double norm(const std::vector<double>& x);
	static std::vector<std::vector<double> > eye(const unsigned int& n);
	static std::vector<std::vector<double> > outerProduct(const std::vector<double>& x,
		const std::vector<double>& y);
};



#endif