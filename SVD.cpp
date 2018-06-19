#include "SVD.hpp"
#include <algorithm>
void SVD::decomp(
	std::vector<std::vector<double> >&U,
	std::vector<std::vector<double> >&S,
	std::vector<std::vector<double> >&V,
	const std::vector<std::vector<double> >&A
)
{
	const double tol = 1e-11;
	unsigned int n = 0, m = 0;

	if (!MXOP::size(n, m, A))
	{
		std::cout << "Invalid input matrix" << std::endl;
		return;
	}

	const int maxIter = 100 * std::max(m, n);


	double err = std::numeric_limits<double>::max();
	int loopCounter = 0;
	U = MXOP::eye(n);
	V = MXOP::eye(m);
	S = MXOP::transpose(A);

	while ((err > tol) && (loopCounter < maxIter) )
	{
		std::vector<std::vector<double> > QU, QV;
		QR::decomp(QU,S,MXOP::transpose(S));
		U = MXOP::multiply(U, QU);
		QR::decomp(QV, S, MXOP::transpose(S));
		V = MXOP::multiply(V,QV);

		double e = MXOP::norm(MXOP::triu(S));
		double f = MXOP::norm(MXOP::diag(S));
		if (abs(f) < 1e-12)
			f = 1.0;
		err = e / f;
		loopCounter++;
	}

	std::vector<double> s = MXOP::diag(S);
	std::vector<std::vector<double> > fixedS = MXOP::zeros(n, m);
	for (unsigned int i = 0; i < s.size(); i++)
	{
		fixedS[i][i] = abs(s[i]);
		if (s[i] < 0)
			for (unsigned int j = 0; j < n; j++)
				U[j][i] *= -1;
	}
	S = fixedS;
}