#include "QR.hpp"
#include "SVD.hpp"
#include <iostream>

int main()
{
	std::vector<std::vector<double> > testMatrixA =
		{
			{ 16,	2,	3,	13	},
			{ 5,	11,	10,	8	},
			{ 9,	7,	6,	12	}
		};

	std::vector<std::vector<double> > testMatrixB =
		{
			{ 17,	24,	1,	8,	15	},
			{ 3,	5,  7,	14,	16	},
			{ 4,	6,	13,	20,	22	},
			{ 10,	12, 19,	21,	3	},
			{ 11,	18,	25,	2,	9	}
		};

	std::vector<std::vector<double> > testMatrixC =
		{
			{ 17,	24,	1	},
			{ 3,	5,  7	},
			{ 4,	6,	13	},
			{ 10,	12, 19	},
			{ 11,	18,	25	}
		};

	std::vector< std::vector< std::vector<double> > > allMatrices =
		{
			testMatrixA, testMatrixB, testMatrixC
		};

	for(const auto& matrix : allMatrices)
	{

		// Result matrices for QR
		std::vector<std::vector<double> > R;
		std::vector<std::vector<double> > Q;

		// Result matrices for SVD
		std::vector<std::vector<double> > U;
		std::vector<std::vector<double> > S;
		std::vector<std::vector<double> > V;

		std::cout << "--------------------------------------" << std::endl;
		std::cout << "QR DECOMPOSITION" << std::endl;
		QR::decomp(Q, R, matrix);

		// Print outputs
		std::cout << "Original matrix" << std::endl;
		MXOP::printMatrix(matrix);
		std::cout << std::endl << "Q: " << std::endl;
		MXOP::printMatrix(Q);
		std::cout << std::endl << "R: " << std::endl;
		MXOP::printMatrix(R);
		std::cout << std::endl << "Q * Q': " << std::endl;
		MXOP::printMatrix(MXOP::multiply(Q, MXOP::transpose(Q)));
		std::cout << std::endl << "Q * R" << std::endl;
		MXOP::printMatrix(MXOP::multiply(Q, R));

		std::cout << "--------------------------------------" << std::endl;
		std::cout << "SVD DECOMPOSITION" << std::endl;
		SVD::decomp(U,S,V,matrix);
		std::cout << "Original matrix" << std::endl;
		MXOP::printMatrix(matrix);
		std::cout << std::endl << "U: " << std::endl;
		MXOP::printMatrix(U);
		std::cout << std::endl << "S: " << std::endl;
		MXOP::printMatrix(S);
		std::cout << std::endl << "V: " << std::endl;
		MXOP::printMatrix(V);
		std::cout << std::endl << "U * S * V': " << std::endl;
		MXOP::printMatrix(MXOP::multiply(U, MXOP::multiply(S,MXOP::transpose(V))));
		std::cout << std::endl << "V * V': " << std::endl;
		MXOP::printMatrix(MXOP::multiply(V, MXOP::transpose(V)));
		std::cout << std::endl << "U * U': " << std::endl;
		MXOP::printMatrix(MXOP::multiply(U, MXOP::transpose(U)));

		std::cout << "--------------------------------------" << std::endl;
		std::cout << "SVD ONLY U" << std::endl;
		SVD::decomp(U, S, V, matrix, SVD::SvdOptions::ONLY_U);
		std::cout << "Original matrix" << std::endl;
		MXOP::printMatrix(matrix);
		std::cout << std::endl << "U: " << std::endl;
		MXOP::printMatrix(U);
		std::cout << std::endl << "S: " << std::endl;
		MXOP::printMatrix(S);
		std::cout << std::endl << "U * U': " << std::endl;
		MXOP::printMatrix(MXOP::multiply(U, MXOP::transpose(U)));

		std::cout << "--------------------------------------" << std::endl;
		std::cout << "SVD ONLY V" << std::endl;
		SVD::decomp(U, S, V, matrix,SVD::SvdOptions::ONLY_V);
		std::cout << "Original matrix" << std::endl;
		MXOP::printMatrix(matrix);
		std::cout << std::endl << "V: " << std::endl;
		MXOP::printMatrix(V);
		std::cout << std::endl << "V * V': " << std::endl;
		MXOP::printMatrix(MXOP::multiply(V, MXOP::transpose(V)));
	}

	return 0;
}