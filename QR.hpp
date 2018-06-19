/********************************************************************
 ********************************************************************
 ** C++ class implementation of a basic QR Decomposition by Axel, 
 ** 2018. 
 **
 **
 ** Performs QR-decomposition of a of an input matrix A, such that
 ** the product of the output matrices, Q and R equals A.
 ** Q is orthonormal, that is, transpose(Q) = inverse(Q),
 ** furhtermore, R, is upper triangular.
 ** 
 **	Example: See Example.cpp
 **
 ** ** This file may be freely copied and distributed! **
 **
 **
 ** This file is distributed in the hope that it will be useful,
 ** but WITHOUT ANY WARRANTY; without even the implied 
 ** warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 ** PURPOSE.  
 **
 ********************************************************************
 ********************************************************************/
#ifndef QR_H
#define QR_H
#include <vector>
#include <math.h>
#include <iostream>
#define QR_SIGN(x) (x < 0) ? -1 : 1
using namespace std;

class QR
{

public:
    
    enum QRSTATUS
    {
        QR_SUCCES,
        QR_ERROR
    };

    static QRSTATUS decomp(
		vector<vector<double> >& outQ,
		vector<vector<double> >& outR,
        const vector<vector<double> >& inA
    );
	static vector<vector<double> > zeros(const unsigned int& n, const unsigned int& m);
	static void printMatrix(const vector<vector<double> >& matrix);
	static void printVector(const vector<double>& vector);
	static vector<vector<double> > multiply(const vector<vector<double> >& A, const vector<vector<double> >& B);
	static vector<vector<double> > transpose(const vector<vector<double> >& mat);
	static double norm(const vector<double>& x);
	static vector<vector<double> > eye(const unsigned int& n);

private:
    static vector<vector<double> > outerProduct(const vector<double>& x,
     const vector<double>& y);
};
#endif