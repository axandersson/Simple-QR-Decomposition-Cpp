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
#include "MatrixOperations.hpp"
#define QR_SIGN(x) ((x < 0) ? -1 : 1)

class QR
{

public:

	enum QRSTATUS
	{
		QR_SUCCES,
		QR_ERROR
	};

	static QRSTATUS decomp(
		std::vector<std::vector<double> >& outQ,
		std::vector<std::vector<double> >& outR,
		const std::vector<std::vector<double> >& inA
	);


};
#endif