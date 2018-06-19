#ifndef SVD_H
#define SVS_H
#include <vector>
#include "MatrixOperations.hpp"
#include "QR.hpp"
class SVD
{
public:
	static void decomp(
		std::vector<std::vector<double> >&outU,
		std::vector<std::vector<double> >&outS,
		std::vector<std::vector<double> >&outV,
		const std::vector<std::vector<double> >&A
	);
};
#endif // !SVD_H
