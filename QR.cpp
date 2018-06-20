#include "QR.hpp"

QR::QRSTATUS QR::decomp(
	std::vector<std::vector<double> >& outQ,
	std::vector<std::vector<double> >& outR,
    const std::vector<std::vector<double> >& inA,
	const int& onlyR
)
{
	// Get matrix size
	unsigned int m, n;
	if (!MXOP::size(m, n, inA))
	{
		return QR::QRSTATUS::QR_ERROR;
	}

	// Number of iterations required for upper-diagonal
    const unsigned int nLoop = std::min(m,n);

	// Initialize output
    outR = inA;
	if (!onlyR)
	{
		outQ = MXOP::eye(m);
	}
   
    for(unsigned int k = 0; k < nLoop; k++)
    {
        std::vector<double> x(m,0.0);
        for(unsigned int j = k; j < m; j++)
        {
            x[j] = outR[j][k];
        }

		x[k] += (QR_SIGN(x[k])) * MXOP::norm(x);
		
		const double normx2 = MXOP::norm(x);
		for(auto& xi : x)
            xi /= normx2;
		
		if (!onlyR)
		{
			outQ = multiplyHousehold(outQ, x);
		}
		outR = multiplyHousehold(outR, x);
    }
    outQ = MXOP::transpose(outQ);
	return QR::QRSTATUS::QR_SUCCES;
}

std::vector<std::vector<double> > QR::multiplyHousehold(std::vector<std::vector<double> > A, const std::vector<double>& x)
{
	 std::vector<double> xa(A[0].size(), 0.0);
	 // -2 * row x, multiplied with A
	 for (unsigned int i = 0; i < xa.size(); i++)
	 {
		 for (unsigned int j = 0; j < A.size(); j++)
		 {
			 xa[i] += x[j]*A[j][i];
		 }
		 xa[i] *= -2.0;
	 }

	 std::vector<std::vector<double> > result = MXOP::outerProduct(x,xa);

	 for (unsigned int i = 0; i < A.size(); i++)
	 {
		 for (unsigned int j = 0; j < A[0].size(); j++)
		 {
			 result[i][j] += A[i][j];
			 if (abs(result[i][j]) < 1e-14)
				 result[i][j] = 0;
		 }
	 }
	 return result;
}
