#include "QR.hpp"

QR::QRSTATUS QR::decomp(
	std::vector<std::vector<double> >& outQ,
	std::vector<std::vector<double> >& outR,
    const std::vector<std::vector<double> >& inA
)
{
    if(inA.empty())
    {
        return QR::QRSTATUS::QR_ERROR;
    }

    if(inA[0].empty())
    {
        return QR::QRSTATUS::QR_ERROR;
    }

    const unsigned int m = inA.size();
    const unsigned int n = inA[0].size();
    const unsigned int nLoop = n < m ? n : m;
    outR = inA;
    outQ = MXOP::eye(m);
    const std::vector<std::vector<double> > identity = MXOP::eye(m);
    for(unsigned int k = 0; k < nLoop; k++)
    {

        std::vector<double> x(m,0.0);
        for(unsigned int j = k; j < m; j++)
        {
            x[j] = outR[j][k];
        }
        double normx = MXOP::norm(x);
        x[k] += (QR_SIGN(x[k])) * normx;
		double normx2 = MXOP::norm(x);
		for(auto& xi : x)
            xi /= normx2;
		std::vector<std::vector<double> > U = MXOP::eye(m);
		std::vector<std::vector<double> > X = MXOP::outerProduct(x,x);
        for (unsigned int i = 0; i < m ; i++)
        {
            for( unsigned int j = 0; j < m; j++)
            {
                U[i][j] -= 2.0 * X[i][j];
            }
        }
        outQ = MXOP::multiply(U,outQ);
        outR = MXOP::multiply(U,outR);
    }
    outQ = MXOP::transpose(outQ);
	return QR::QRSTATUS::QR_SUCCES;
}
