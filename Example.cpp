#include "QR.hpp"
#include <iostream>

int main()
{
	// Matrix to be decomposed
	vector<vector<double> > matrix = 
	{	
		{ 16,2,3,13},
		{ 5,11,10,8},
		{ 9,7,6,12}
	};

	// Result matrices
	vector<vector<double> > R;
	vector<vector<double> > Q;

	QR::decomp(Q, R, matrix);

	// Print outputs
	cout << "Original matrix" << endl;
	QR::printMatrix(matrix);
	cout << "Q: " << endl;
	QR::printMatrix(Q);
	cout << "R: " << endl;
	QR::printMatrix(R);
	cout << "Q * Q': " << endl;
	QR::printMatrix(QR::multiply(Q, QR::transpose(Q)));
	cout << "Q * R" << endl;
	QR::printMatrix(QR::multiply(Q, R));
	return 0;
}