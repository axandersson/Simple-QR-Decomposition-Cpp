#include "MatrixOperations.hpp"
#include "QR.hpp"

std::vector<double> MXOP::triu(const std::vector<std::vector<double> >& A)
{
	const unsigned int n = A.size();
	const unsigned int m = A[0].size();
	std::vector<double> upperDiagonal;
	const unsigned int minMN= std::min(n,m);
	for (unsigned int i = 0; i < minMN; i++)
		for(unsigned int j=i+1;j<m;j++)
			upperDiagonal.push_back(A[i][j]);
	return upperDiagonal;
}


std::vector<double> MXOP::diag(const std::vector<std::vector<double> >& A)
{
	std::vector<double> diagonal;
	const unsigned int n = std::min(A.size(), A[0].size());
	for (unsigned int i = 0; i < n; i++)
		diagonal.push_back(A[i][i]);
	return diagonal;
}

bool MXOP::size(unsigned int& n, unsigned int& m, const std::vector<std::vector<double> >& matrix)
{
	n = matrix.size();
	if (n == 0)
		return false;
	m = matrix[0].size();
	if (m == 0)
		return false;
	return true;
}

void MXOP::printMatrix(const std::vector<std::vector<double> >& matrix)
{
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		for (unsigned int j = 0; j < matrix[i].size(); j++)
		{
			std::cout << matrix[i][j] << '\t';
		}
		std::cout << std::endl;
	}
}

void MXOP::printVector(const std::vector<double>& vec)
{
	for (unsigned int j = 0; j < vec.size(); j++)
	{
		std::cout << vec[j] << ' ';
	}
	std::cout << std::endl;

}
std::vector<std::vector<double> > MXOP::zeros(
	const unsigned int& n,
	const unsigned int& m
)
{
	std::vector<std::vector<double> > vec(n, std::vector<double>(m, 0.0));
	return vec;
}
std::vector<std::vector<double> > MXOP::transpose(const std::vector<std::vector<double> >& mat)
{
	const unsigned int n = mat.size();
	const unsigned int m = mat[0].size();
	std::vector<std::vector<double> > trans = zeros(m, n);
	for (unsigned int i = 0; i < n; i++)
	{
		for (unsigned int j = 0; j < m; j++)
		{
			trans[j][i] = mat[i][j];
		}
	}
	return trans;
}

double MXOP::norm(const std::vector<double>& x)
{
	double n = 0.0;
	for (const auto& xi : x)
		n += xi * xi;
	return sqrt(n);
}
std::vector<std::vector<double> > MXOP::multiply(
	const std::vector<std::vector<double> >& A,
	const std::vector<std::vector<double> >&B
)
{
	const unsigned int nA = A.size();
	const unsigned int mA = A[0].size();
	const unsigned int mB = B[0].size();
	std::vector<std::vector<double> > out = zeros(nA, mB);
	for (unsigned int i = 0; i < nA; i++)
	{
		for (unsigned int j = 0; j < mB; j++)
		{
			for (unsigned int k = 0; k < mA; k++)
			{
				out[i][j] += A[i][k] * B[k][j];
			}
			if (abs(out[i][j]) < 1e-14)
			{
				out[i][j] = 0;
			}
		}
	}
	return out;
}


std::vector<std::vector<double> > MXOP::outerProduct(
	const std::vector<double>& x,
	const std::vector<double>& y
)
{
	std::vector<std::vector<double> > outer = zeros(x.size(), y.size());
	for (unsigned int i = 0; i < x.size(); i++)
	{
		for (unsigned int j = 0; j < y.size(); j++)
		{
			outer[i][j] = x[i] * y[j];
		}
	}
	return outer;
}
std::vector<std::vector<double> > MXOP::eye(const unsigned int& n)
{
	std::vector<std::vector<double> > ident = zeros(n, n);
	for (unsigned int i = 0; i <n; i++)
		ident[i][i] = 1.0;
	return ident;
}
