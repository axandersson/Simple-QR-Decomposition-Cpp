#include "QR.hpp"

QR::QRSTATUS QR::decomp(
	vector<vector<double> >& outQ,
	vector<vector<double> >& outR,
    const vector<vector<double> >& inA
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
    outQ = eye(m);
    const vector<vector<double> > identity = eye(m);
    for(unsigned int k = 0; k < nLoop; k++)
    {

        vector<double> x(m,0.0);
        for(unsigned int j = k; j < m; j++)
        {
            x[j] = outR[j][k];
        }
        double normx = norm(x);
        x[k] += (QR_SIGN(x[k])) * normx;
		double normx2 = norm(x);
		for(auto& xi : x)
            xi /= normx2;
		vector<vector<double> > U = eye(m);
		vector<vector<double> > X = outerProduct(x,x);
        for (unsigned int i = 0; i < m ; i++)
        {
            for( unsigned int j = 0; j < m; j++)
            {
                U[i][j] -= 2.0 * X[i][j];
            }
        }
        outQ = multiply(U,outQ);
        outR = multiply(U,outR);
    }
    outQ = transpose(outQ);
	return QR::QRSTATUS::QR_SUCCES;
}

void QR::printMatrix(const vector<vector<double> >& matrix)
{
	for (unsigned int i = 0; i < matrix.size(); i++)
	{
		for (unsigned int j = 0; j < matrix[i].size(); j++)
		{
			cout << matrix[i][j] << ' ';
		}
		cout << endl;
	}
}

void QR::printVector(const vector<double>& vec)
{
	for (unsigned int j = 0; j < vec.size(); j++)
	{
		cout << vec[j] << ' ';
	}
	cout << endl;

}
vector<vector<double> > QR::zeros(
	const unsigned int& n, 
	const unsigned int& m
)
{
	vector<vector<double> > vec(n, vector<double>(m, 0.0));
	return vec;
}
vector<vector<double> > QR::transpose(const vector<vector<double> >& mat)
{
    const unsigned int n = mat.size();
    const unsigned int m = mat[0].size();
	vector<vector<double> > trans = zeros(m,n);
    for(unsigned int i = 0; i < n; i++)
    {
        for(unsigned int j = 0; j < m; j++)
        {
            trans[j][i] = mat[i][j];
        }
    }
    return trans;
}

double QR::norm(const vector<double>& x)
{
    double n = 0.0;
    for(const auto& xi : x)
        n += xi*xi;
    return sqrt(n);
}
vector<vector<double> > QR::multiply(
	const vector<vector<double> >& A, 
	const vector<vector<double> >&B
)
{
    const unsigned int nA = A.size();
    const unsigned int mA = A[0].size();
    const unsigned int mB = B[0].size();
	vector<vector<double> > out = zeros(nA,mB);
    for(unsigned int i = 0; i < nA; i++)
    {
        for(unsigned int j = 0; j < mB; j++)
        {
            for(unsigned int k = 0; k < mA; k++)
            {
				double prod = A[i][k] * B[k][j];
				if (abs(prod) < 1e-11)
				{
					prod = 0;
				}
				out[i][j] += prod;
            }
			if (abs(out[i][j]) < 1e-11)
				out[i][j] = 0;
        }
    }
	return out;
}


vector<vector<double> > QR::outerProduct(
	const vector<double>& x,
	const vector<double>& y
)
{
	vector<vector<double> > outer = zeros(x.size(),y.size());
    for(unsigned int i = 0; i < x.size(); i++)
    {
        for(unsigned int j=0; j < y.size();j++) 
        {
            outer[i][j] = x[i]*y[j];
        }
    }
	return outer;
}
vector<vector<double> > QR::eye(const unsigned int& n)
{
	vector<vector<double> > ident = zeros(n,n);
    for(unsigned int i = 0; i <n; i++)
        ident[i][i] = 1.0;
    return ident;
}