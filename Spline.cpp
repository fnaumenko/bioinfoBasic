#include "Spline.h"


/************************ Bezier2D ************************/

// Bezier 2D curve
// https://www.codeproject.com/Articles/25237/Bezier-Curves-Made-Simple
static class Bezier2D
{
	static const float factorials[];	// factorials 'table'

	static float factorial(int n) {
		assert(n >= 0);
		assert(n <= 32);
		//if (n < 0) { throw range_error(to_string(n) + " is less than 32"); }
		//if (n > 32) { throw range_error(to_string(n) + " is greater than 32"); }
		return factorials[n];
	}

	// Calculate Bernstein basis
	//	@param n: number of points
	//	@param i: point index
	//	@param t: distance
	static float Bernstein(int n, int i, float t) {
		return factorial(n) / (factorial(i) * factorial(n - i))
			* (!i && !t ? 1.0 : pow(t, i))						// t^i
			* (n == i && t == 1 ? 1.0 : pow((1 - t), (n - i)));	// (1 - t)^i
	}

public:
	static void RefineSummit(coviter& it0, int32_t halfBase, point& summit
#ifdef DEBUG_OUTPUT
		, const char* chrName, ofstream& ostream
#endif
	) {
		assert(halfBase < 16);
		float val = 0;

		// Calculate splined curve
		const auto itCntDecr = 2 * halfBase;	// number of iterators within splined range, decremented by 1
		const auto itCnt = itCntDecr + 1;		// number of iterators within splined range

		advance(it0, halfBase);			// now it0 points to the end of splined range!
		chrlen ptCnt = it0->first;		// last pos in range - temporary
		advance(it0, -2 * halfBase);	// now it0 points to the start of splined range!
		ptCnt -= it0->first;			// number of points within splined range

		const float	step = 1.0 / (ptCnt - 1);
		const auto posEnd = it0->first + ptCnt;
		float t = 0;
		summit.second = 0;
		// loop loop over points in range
		for (auto pos = it0->first; pos < posEnd;) {
			auto it = it0;				// start from the beginning of the region
			val = 0;
			// loop over iterators in range
			for (int8_t i = 0; i < itCnt; i++, it++)
				val += it->second * Bernstein(itCntDecr, i, t);	// splined value
			t += step;
			if ((1.0 - t) < 5e-6)	t = 1.0;

			// find refined summit
			if (val > summit.second)
				summit.first = pos,
				summit.second = val;
#ifdef DEBUG_OUTPUT
			if (ostream.is_open()) {
				ostream << chrName << TAB << pos;
				ostream << TAB << ++pos << TAB << val << LF;
			}
			else
#endif
				pos++;
		}
	}

} bezier2D;

// factorials 'table'
const float Bezier2D::factorials[] = {
	1.0,									// 0
	1.0,									// 1
	2.0,									// 2
	6.0,									// 3
	24.0,									// 4
	120.0,									// 5
	720.0,									// 6
	5040.0,									// 7 
	40320.0,								// 8 
	362880.0,								// 9 
	3628800.0,								// 10 
	39916800.0,								// 11
	479001600.0,							// 12
	6227020800.0,							// 13
	87178291200.0,							// 14
	1307674368000.0,						// 15
	20922789888000.0,						// 16
	355687428096000.0,						// 17
	6402373705728000.0,						// 18
	121645100408832000.0,					// 19
	2432902008176640000.0,					// 20
	51090942171709440000.0,					// 21
	1124000727777607680000.0,				// 22
	25852016738884976640000.0,				// 23
	620448401733239439360000.0,				// 24
	15511210043330985984000000.0,			// 25
	403291461126605635584000000.0,			// 26
	10888869450418352160768000000.0,		// 27
	304888344611713860501504000000.0,		// 28
	8841761993739701954543616000000.0,		// 29
	265252859812191058636308480000000.0,	// 30
	8222838654177922817725562880000000.0,	// 31
	263130836933693530167218012160000000.0	// 32
};

/************************ Bezier2D: end ************************/

