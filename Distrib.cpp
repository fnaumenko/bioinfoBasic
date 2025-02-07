/**********************************************************
Distrib.cpp
Last modified: 02/07/2025
***********************************************************/

#include "Distrib.h"
#include "spline.h"
#include "DataReader.h"
#include <algorithm>    // std::sort

// square of doubled Pi
const float SDPI = float(sqrt(3.1415926 * 2));
// ratio of the summit height to height of the measuring point
const int hRatio = 2;
// log of ratio of the summit height to height of the measuring point
const float lghRatio = float(log(hRatio));

// Approximate distribution formulas
struct ADF {
	const char* Title;

	// Returns two common factors of the distrib equation
	//	@param p: distrib params: mean/alpha and sigma/beta
	fpair (*CommonFactors)(const fpair& p);

	// Returns the value of the distribution function at a given point and common factors
	//	@param p: distrib params: mean/alpha and sigma/beta
	//	@param x: x-coordinate
	//	@param commFactors: two common factors of the distrib equation
	double (*Value)	(const fpair& p, fraglen x, const fpair& commFactors);

	float (*Mode)	(const fpair& p);
	float (*Mean)	(const fpair& p);
	float (*Median)	(const fpair& p);

	// Calculates approximate distribution parameters
	//	@param keypts[in]: key points: X-coord of highest point, X-coord of half height of the highest point
	//	@param p[out]: returned params: mean(alpha) & sigma(beta)
	void (*CalcParams)(const fpair& keypts, fpair& p);

	// Returns the value of the distribution function at a given point
	//	@param p: distrib params: mean/alpha and sigma/beta
	//	@param x: x-coordinate
	double GetValue	(const fpair& p, fraglen x) { return Value(p, x, CommonFactors(p)); }
};

static ADF ADFs[Distrib::eDType::CNT] {
#define FACTOR1	cmFactors.first
#define FACTOR2	cmFactors.second
#define MEAN	p.first
#define SIGMA	p.second
#define ALPHA	p.first
#define BETA	p.second
#define POW2_SIGMA	SIGMA * SIGMA
#define X_HIGH		keypts.first
#define X_HALFHIGH	keypts.second

	{ "Norm",
	[](const fpair& p) ->fpair { return { SIGMA * SDPI, 0.f}; },	// CommonFactors
	[](const fpair& p, fraglen x, const fpair& cmFactors) { // Value
		return exp(-pow(((x - MEAN) / SIGMA), 2) / 2) / FACTOR1;
	},
	[](const fpair& p) { return MEAN; },					// Mode
	[](const fpair& p) { return MEAN; },					// Mean
	[](const fpair& p) { return MEAN; },					// Median
	[](const fpair& keypts, fpair& p) {						// CalcParams
		MEAN = X_HIGH;
		SIGMA = float(sqrt(pow(X_HALFHIGH - MEAN, 2) / lghRatio / 2));
	},
	},
	{ "Lognorm",
	[](const fpair& p) ->fpair { return { SIGMA * SDPI, 2 * POW2_SIGMA}; },	// CommonFactors
	[](const fpair& p, fraglen x, const fpair& cmFactors) { 	// Value
		return exp(-pow((log(x) - MEAN), 2) / FACTOR2) / (FACTOR1 * x);
	},
	[](const fpair& p) { return exp(MEAN - POW2_SIGMA); },		// Mode
	[](const fpair& p) { return exp(MEAN + POW2_SIGMA / 2); },	// Mean
	[](const fpair& p) { return exp(MEAN); },					// Median
	[](const fpair& keypts, fpair& p) {							// CalcParams
		const float lgM = log(X_HIGH);		// logarifm of Mode
		const float lgH = log(X_HALFHIGH);	// logarifm of middle height
		MEAN = (lgM * (lghRatio + lgM - lgH) + (lgH * lgH - lgM * lgM) / 2) / lghRatio;
		SIGMA = sqrt(MEAN - lgM);
	},
	},
	{ "Gamma",
	[](const fpair& p) ->fpair { return { ALPHA - 1, float(pow(BETA, ALPHA)) }; },	// CommonFactors
	[](const fpair& p, fraglen x, const fpair& cmFactors) { // Value
		return pow(x, FACTOR1) * exp(-(x / BETA)) / FACTOR2;
	},
	[](const fpair& p) { return (ALPHA - 1) * BETA; },		// Mode
	[](const fpair& p) { return ALPHA * BETA; },			// Mean 
	[](const fpair& p) { return 0.f; },						// Median: undefined
	[](const fpair& keypts, fpair& p) {						// CalcParams
		BETA = (X_HALFHIGH - X_HIGH * (1 + log(X_HALFHIGH / X_HIGH))) / lghRatio;
		ALPHA = (X_HIGH / BETA) + 1;
	}
	},
};

// Sets half height X-coordinate on the right slope
//	@param p0[in,out]: previous point
//	@param p[in,out]: current point
//	@param summit[in,out]: summit point
//	@param halfSummitX[in,out]: returned value
void SetHalfSummitX(fpair& p0, fpair& p, fpair& summit, float& halfSummitX)
{
#ifdef PRINT
	std::printf("%.2f\t%.2f\n", p.first, p.second);
#endif
	if (p.second >= summit.second)
		p.swap(summit);
	else {
		if (p.second < summit.second / hRatio) {
			if (!halfSummitX)
				halfSummitX = p0.first + p0.second / (p.second + p0.second);
		}
		p.swap(p0);
	}
}

//===== Bezier2D
// 
// Bezier 2D curve
// https://www.codeproject.com/Articles/25237/Bezier-Curves-Made-Simple
static class Bezier2D
{
public:
	static const BYTE MAX_POINT_CNT = 98;

	// Performes Bezier interpolation and return the position of the maximum of the Bezier curve
	//	@param pts: raw points
	//	@param outPtCnt: number of Bezier curve points
	//	@returns: position of the maximum of the Bezier curve
	static fpair GetKeyPoints(const dmap& pts, float cutoffThreshold)
	{
		map<int, chrlen>::const_iterator it0;	// start it
		map<int, chrlen>::const_iterator it1;	// end it
		const auto ptCnt = Trim(pts, it0, it1, cutoffThreshold);

		if (ptCnt > MAX_POINT_CNT)
			throw range_error("Bezier2D: number of points " + to_string(ptCnt) + " is greater than maximum permissible " + to_string(MAX_POINT_CNT));
		const USHORT outPtCnt = it1->first - it0->first;
		const float	step = 1.f / (outPtCnt - 1);
		float	d = 0;								// distance
		fpair summit{};
		fpair p0;
		float halfSummitX = 0;

#ifdef PRINT
		printf("count: %u  skipCnt: %zu\n", ptCnt + 1, pts.size() - ptCnt - 1);
		printf("BEZIER SPLINED DIFFS FREQUENCY  %d\n", outPtCnt);
#endif
		++it1;
		// Calculate points on curve
		for (UINT pInd = 0; pInd < outPtCnt; pInd++) {
			if ((1.f - d) < 5e-6)
				d = 1.f;
			fpair p;		// interpolated point
			BYTE i = 0;
			for (auto it = it0; it != it1; it++) {
				auto basis = Bernstein(ptCnt, i++, d);
				p.first += basis * it->first;
				p.second += basis * it->second;
			}
			d += step;

			SetHalfSummitX(p0, p, summit, halfSummitX);
#ifndef PRINT
			if (halfSummitX)
				break;
#endif
		}
#ifdef PRINT
		printf("MAX POS: %.1f  HALF POS %.1f:\n", summit.first, halfSummitX);
#endif
		return fpair(
			summit.first,							// summit X-coord
			halfSummitX
		);
	}

private:
	// factorials 'table'
	static const double factorials[MAX_POINT_CNT + 1];

	static BYTE Trim(const dmap& pts,
		map<int, chrlen>::const_iterator& it0,
		map<int, chrlen>::const_iterator& it1,
		float cutoffThreshold)
	{
		it0 = pts.begin();		// start it
		it1 = prev(pts.end());	// end it
		if (pts.size() <= 5)	return BYTE(pts.size());

		// ** cut off single frequency iterators at the edges
		// ** trim the distribution's 'tails'
		// define max value
		float maxVal = 0;
		for (const auto& f : pts)
			if (maxVal < f.second)
				maxVal = f.second;

		// trim entries with value less than cutoffThreshold of max value
		maxVal *= cutoffThreshold;
		while (it0->second < maxVal)	it0++;
		it0--;
		while (it1->second < maxVal)	it1--;
		it1++;
		return BYTE(distance(it0, it1));
	}

	// Calculate Bernstein basis
	//	@param ptCnt: number of points
	//	@param ptInd: point index
	//	@param d: distance
	//	@returns: Bernstein basis
	static float Bernstein(BYTE ptCnt, BYTE ptInd, float d)
	{
		// Prevent problems with pow
		float ti = !d && !ptInd ? 1.f : float(pow(d, ptInd));	// d^i
		float xi = ptCnt == ptInd && d == 1.f ?					// (1 - d)^i
			1.f :
			float(pow((1 - d), (ptCnt - ptInd)));

		double a1 = factorials[ptCnt];
		double a2 = factorials[ptInd];
		double a3 = factorials[ptCnt - ptInd];

		return ti * xi * float(a1 / (a2 * a3));
	}

} bezier2D;

const double Bezier2D::factorials[] = {
	1.,	// 0!
	1.,
	2.,	// 2!
	6.,
	24.,	// 4!
	120.,
	720.,	// 6!
	5040.,
	40320.,	// 8!
	362880.,
	3628800.,	// 10!
	39916800.,
	479001600.,	// 12!
	6227020800.,
	87178291200.,	// 14!
	1307674368000.,
	20922789888000.,	// 16!
	355687428096000.,
	6402373705728000.,	// 18!
	121645100408832000.,
	2432902008176640000.,	// 20!
	51090942171709440000.,
	1124000727777607680000.,	// 22!
	25852016738884976640000.,
	620448401733239439360000.,	// 24!
	15511210043330985984000000.,
	403291461126605635584000000.,	// 26!
	10888869450418352160768000000.,
	304888344611713860501504000000.,	// 28!
	8841761993739701954543616000000.,
	2.6525285981219105863630848e+32,	// 30!
	8.22283865417792281772556288e+33,
	2.6313083693369353016721801216e+35,	// 32!
	8.68331761881188649551819440128e+36,
	2.9523279903960414084761860964352e+38,	// 34!
	1.0333147966386144929666651337523e+40,
	3.7199332678990121746799944815084e+41,	// 36!
	1.3763753091226345046315979581581e+43,
	5.2302261746660111176000722410007e+44,	// 38!
	2.0397882081197443358640281739903e+46,
	8.1591528324789773434561126959612e+47,	// 40!
	3.3452526613163807108170062053441e+49,
	1.4050061177528798985431426062445e+51,	// 42!
	6.0415263063373835637355132068514e+52,
	2.6582715747884487680436258110146e+54,	// 44!
	1.1962222086548019456196316149566e+56,
	5.5026221598120889498503054288003e+57,	// 46!
	2.5862324151116818064296435515361e+59,
	1.2413915592536072670862289047373e+61,	// 48!
	6.082818640342675608722521633213e+62,
	3.0414093201713378043612608166065e+64,	// 50!
	1.5511187532873822802242430164693e+66,
	8.0658175170943878571660636856404e+67,	// 52!
	4.2748832840600255642980137533894e+69,
	2.3084369733924138047209274268303e+71,	// 54!
	1.2696403353658275925965100847567e+73,
	7.1099858780486345185404564746372e+74,	// 56!
	4.0526919504877216755680601905432e+76,
	2.3505613312828785718294749105151e+78,	// 58!
	1.3868311854568983573793901972039e+80,
	8.3209871127413901442763411832234e+81,	// 60!
	5.0758021387722479880085681217663e+83,
	3.1469973260387937525653122354951e+85,	// 62!
	1.9826083154044400641161467083619e+87,
	1.2688693218588416410343338933516e+89,	// 64!
	8.2476505920824706667231703067855e+90,
	5.4434493907744306400372924024784e+92,	// 66!
	3.6471110918188685288249859096605e+94,
	2.4800355424368305996009904185692e+96,	// 68!
	1.7112245242814131137246833888127e+98,
	1.1978571669969891796072783721689e+100,	// 70!
	8.5047858856786231752116764423993e+101,
	6.1234458376886086861524070385275e+103,	// 72!
	4.4701154615126843408912571381251e+105,
	3.3078854415193864122595302822125e+107,	// 74!
	2.4809140811395398091946477116594e+109,
	1.8854947016660502549879322608611e+111,	// 76!
	1.4518309202828586963407078408631e+113,
	1.1324281178206297831457521158732e+115,	// 78!
	8.9461821307829752868514417153983e+116,
	7.1569457046263802294811533723187e+118,	// 80!
	5.7971260207473679858797342315781e+120,
	4.753643337012841748421382069894e+122,	// 82!
	3.9455239697206586511897471180121e+124,
	3.3142401345653532669993875791301e+126,	// 84!
	2.8171041143805502769494794422606e+128,
	2.4227095383672732381765523203441e+130,	// 86!
	2.1077572983795277172136005186994e+132,
	1.8548264225739843911479684564555e+134,	// 88!
	1.6507955160908461081216919262454e+136,
	1.4857159644817614973095227336208e+138,	// 90!
	1.352001527678402962551665687595e+140,
	1.2438414054641307255475324325874e+142,	// 92!
	1.1567725070816415747592051623062e+144,
	1.0873661566567430802736528525679e+146,	// 94!
	1.0329978488239059262599702099395e+148,
	9.9167793487094968920957140154189e+149,	// 96!
	9.6192759682482119853328425949564e+151,
	9.4268904488832477456261857430572e+153,	// 98! : max possible value
	// (!98)^2 = 8.8866263535246200177702174899954e+307, while max double value is 1.7976931348623158e+308
};


//===== Distrib::ADPs

fpair Distrib::ADPs::ADP::GetSplineKeyPoints(const dmap& distr, fraglen base, fpair& summit)
{
	//dpoint p0, p;	// previous, current point
	fpair p0, p;	// previous, current point
	float halfSummitX = 0;
	SSpliner<dVal_t> spliner(
#ifdef MY_DEBUG					// to visualize ROUGH or SMOOTH distributions separately
		eCurveType::ROUGH,
#else
		base <= Distrib::smoothBase ? eCurveType::ROUGH : eCurveType::SMOOTH,
#endif
		base);
#ifdef MY_DEBUG
	_spline.clear();
	fpair keyPts(0, 0);
#endif

	summit.second = 0;
	for (auto& rawp : distr) {
		p.first = spliner.CorrectX(rawp.first);	// X: minus MA & MM base back shift
		p.second = spliner.Push(rawp.second);	// Y: splined

		SetHalfSummitX(p0, p, summit, halfSummitX);
#ifndef PRINT
		if (halfSummitX)
			break;
#endif
	}
	return fpair(
		summit.first,							// summit X-coord
		halfSummitX
	);
}

void Distrib::ADPs::ADP::SetParamsForSpline(Distrib& distr, eDIndex dind, fpair& keypts)
{
	auto calcParams = ADFs[dind].CalcParams;
	const BYTE failCntLim = 2;	// max count of base's decreasing steps after which PCC is considered only decreasing
	BYTE failCnt = 0;	// counter of base's decreasing steps after which PCC is considered only decreasing

	distr.EstimateSplineBase();
#ifdef MY_DEBUG
	int i = 0;					// counter of steps
#endif
	// calculate the highest PCC by iteratively searching through the 'base' values
	for (fraglen base = distr._base; base; base--) {
		//dpoint summit;
		fpair summit;
		ADP adp;
		auto keypts0 = GetSplineKeyPoints(distr, base, summit);

		calcParams(keypts0, adp.Params);
		adp.CalcPCC(distr, dind, keypts0.first);
#ifdef MY_DEBUG
		* _s << setw(4) << setfill(SPACE) << left << ++i;
		*_s << "base: " << setw(2) << base << "  summitX: " << keypts.first << "\tpcc: " << adp.PCC;
		if (adp > adp)	*_s << "\t>";
		*_s << LF;
		if (_fillSpline) { for (dpoint p : _spline)	*_s << p.first << TAB << p.second << LF; _fillSpline = false; }
#endif
		if (adp > *this) {
			std::swap(*this, adp);
			distr._summit.swap(summit);
			keypts.swap(keypts0);
			failCnt = 0;
		}
		else {
			if (adp.PCC > 0)	failCnt++;		// negative PCC is possible in rare cases
			else if (adp.IsUndefPcc()) {
				SetUndefPcc();
				break;
			}
			if (failCnt > failCntLim)
				break;
		}
	}
#ifdef MY_DEBUG
	* _s << LF;
#endif
}

void Distrib::ADPs::ADP::SetParamsForInterpol(Distrib& distr, eDIndex dind, fpair& keypts)
{
	auto calcParams = ADFs[dind].CalcParams;
	ADP adp;

	keypts = Bezier2D::GetKeyPoints(distr, 0.2);
	calcParams(keypts, adp.Params);
	adp.CalcPCC(distr, dind, keypts.first);
}

void Distrib::ADPs::ADP::CalcPCC(const dmap& distr, eDIndex dind, int Mode, bool full)
{
	const fpair commFactors = ADFs[dind].CommonFactors(Params);	// two constant terms of the distrib equation
	const auto dVal = ADFs[dind].Value;	// function that calculates the 'type' distribution coordinate
	const double cutoffY = dVal(Params, Mode, commFactors) / 1000;	// break when Y became less then 0.1% of max value
	double	sumA = 0, sumA2 = 0;	// sum, sum of squares of original values
	double	sumB = 0, sumB2 = 0;	// sum, sum of squares of calculated values
	double	sumAB = 0;				// sum of products of original and calculated values
	UINT	cnt = 0;				// count of points

	// one pass PCC calculation
	SetUndefPcc();
	for (auto& f : distr) {
		if (!full && f.first < Mode)		continue;
		const double b = dVal(Params, f.first, commFactors);	// y-coordinate (value) of the calculated sequence
		if (isNaN(b))						return;
		if (f.first > Mode && b < cutoffY)	break;
		const double a = double(f.second);							// y-coordinate (value) of the original sequence
		sumA += a;
		sumB += b;
		sumA2 += a * a;
		sumB2 += b * b;
		sumAB += a * b;
		cnt++;
	}
	float pcc = float((sumAB * cnt - sumA * sumB) /
		sqrt((sumA2 * cnt - sumA * sumA) * (sumB2 * cnt - sumB * sumB)));
	if (!isNaN(pcc))	PCC = pcc;
}

#define SETW left<<setw(4)
#define UNITAB SETW<<SPACE<<TAB	// tab stretching 4 spaces to display regardless of tab size (4 or 8)

void Distrib::ADPs::IndexedADP::Print(dostream& s, float maxPCC) const
{
	if (IsSet) {
		auto& adfs = ADFs[DIndex];

		s << adfs.Title << TAB;
		if (IsUndefPcc())
			s << "parameters cannot be called";
		else {
			s << setprecision(5) << PCC << TAB;

			// ** print PCC
			if (maxPCC)
				// print percent to max PCC
				if (maxPCC != PCC)
					s << setprecision(3) << 100 * ((PCC - maxPCC) / maxPCC) << "%\t";
				else
					s << UNITAB;

			// ** print basic params
			s << SETW << setprecision(4) << Params.first << TAB << Params.second << TAB;

			// ** print derived params
			s	<< SETW << adfs.Mode(Params) << TAB
				<< SETW << adfs.Mean(Params) << TAB;
			float median = adfs.Median(Params);
			if (median)		// for gamma Median is equal to 0
				s << SETW << median;
		}
		s << LF;
	}
}

bool Distrib::ADPs::IsSetInSorted(eDType dtype) const
{
	const eDIndex dind = GetDType(dtype);
	for (auto& dp : _indADPs)
		if (dp.DIndex == dind)
			return dp.IsSet;
	return false;
}

int Distrib::ADPs::SetSortedCount() const
{
	int cnt = 0;
	for (auto& dp : _indADPs)
		cnt += dp.IsSet;
	return cnt;
}

void Distrib::ADPs::Sort()
{
	if (!_sorted) {
		sort(_indADPs.begin(), _indADPs.end(),
			[](const ADP& dp1, const ADP& dp2) -> bool
			{ return dp1 > dp2; }
		);
		_sorted = true;
	}
}

Distrib::ADPs::ADPs()
{
	int i = 0;
	for (auto& dp : _indADPs)
		dp.DIndex = eDIndex(i++);
}

bool Distrib::ADPs::SetDTypes(eDType dtype)
{
	for (auto& dp : _indADPs)
		dp.IsSet = IsIndex(dtype, dp.DIndex);
	if (_indADPs[1].IsSet && !_indADPs[0].IsSet)
		return _indADPs[0].IsSet = true;
	return false;
}

void Distrib::ADPs::CalcParams(eDType dtype, eDSmooth smode, Distrib& distr)
{
	ADP adp;
	fpair keypts;
	bool checkNormThreshold = SetDTypes(dtype);
	eDIndex dind = _indADPs[iLNORM].IsSet ? iLNORM : (_indADPs[iNORM].IsSet ? iNORM : iGAMMA);

	// set params for the preferred distribution type
	if (smode == eDSmooth::SPLINE)
		adp.SetParamsForSpline(distr, dind, keypts);
	else
		adp.SetParamsForInterpol(distr, dind, keypts);

	CopyParams(dind, adp);
	// set params for other defined distribution types
	for (auto& dp : _indADPs)
		if (dp.IsSet && dp.DIndex != dind) {
			ADFs[dp.DIndex].CalcParams(keypts, adp.Params);
			adp.CalcPCC(distr, dp.DIndex, keypts.first);
			CopyParams(dp.DIndex, adp);
		}

	// clear normal distribution if its PCC is less then lognorm PCC by the threshold 2%
	if (checkNormThreshold && _indADPs[iLNORM].PCC / _indADPs[iNORM].PCC > 1.02F)
		_indADPs[LNORM].IsSet = false;
}

void Distrib::ADPs::Print(dostream& s)
{
	static const char* N[] = { "mean", "sigma" };	// normal, lognormal parameters
	static const char* G[] = { "alpha", "beta" };	// gamma parameters
	static const char* P[] = { "p1", "p2" };		// unified parameters
	static const char* a[] = { "* ", "**" };		// asterisks - footnotes
	const bool notSingle = SetSortedCount() > 1;	// more then 1 output distr type
	float maxPCC = 0;

	Sort();			// should already be sorted by PrintWarning(), but just in case
	const bool isGamma = IsSetInSorted(eDType::GAMMA);

	// ** print title
	s << LF << UNITAB << " PCC\t";
	if (notSingle)
		s << "relPCC\t",
		maxPCC = _indADPs[iNORM].PCC;
	if (!isGamma)		s << N[0] << TAB << N[1];
	else if (notSingle)	s << P[0] << a[0] << TAB << P[1] << a[1];
	else				s << G[0] << TAB << G[1];
	s << "\tMode\tMean";
	if (notSingle || !isGamma)
		s << "\tMedian";

	// ** print values
	s << LF;
	for (const auto& params : _indADPs)
		params.Print(s, maxPCC);

	// ** print note
	if (notSingle && isGamma) {
		auto title = ADFs[2].Title;
		s << LF;
		for (BYTE i = 0; i < 2; i++)
			s << setw(3) << a[i] << P[i] << " - " << N[i] << ", or "
			<< G[i] << " for " << title << LF;
	}
}


//===== Distrib

const char* Distrib::sDistrib = "distribution";

double Distrib::GetApprValue(eDType dtype, float mean, float sigma, fraglen x)
{
	return ADFs[ADPs::GetDType(dtype)].GetValue({mean, sigma}, x);
}

const string Distrib::Spec(eSpec s) { 
	const string sSpec[] = {
		"is degenerate",
		"is smooth",
		"is modulated",
		"is even",
		"is trimmed from the left",
		"is heavily trimmed from the left",
		"looks slightly defective on the left",
		"looks defective on the left"
	};

	return "Distribution " + sSpec[int(s)];
}

const string sParams = "parameters";

void Distrib::EstimateSplineBase()
{
	using rpoint = pair<int, dVal_t>;	// initial raw sequence point
	fraglen halfX = 0;
	rpoint pMin(0, 0), pMMax(pMin);	// minimum, maximum point
	vector<value_type> extr;		// local extremes

	//== define pMin, pMMax and halfX
	{
		const int CutoffFrac = 100;	// fraction of the maximum height below which scanning stops on the first pass
		dVal_t	cutoffY = 0;		// Y-value below which scanning stops on the first pass
		USHORT peakCnt = 0;
		bool up = false;
		auto it = begin();
		SSpliner<dVal_t> spliner(eCurveType::ROUGH, 1);

		_base = 0;
		extr.reserve(20);
		for (rpoint p0(*it++), pMax(pMin); it != end(); it++) {
			rpoint p{
				it->first,
				dVal_t(spliner.Push(it->second))
			};
			if (p.second > p0.second) {		// increasing curve
				if (!up)					// treat pit
					extr.push_back(pMax), up = true, pMin.swap(p0);
			}
			else {							// decreasing curve
				if (up) {					// treat peak
					extr.push_back(pMin), up = false, pMax = p0;

					if (p0.second > pMMax.second) {
						pMMax.swap(p0);
						cutoffY = pMMax.second / CutoffFrac;
					}
					peakCnt++;
				}
				if (peakCnt && p.second >= pMMax.second / 2)
					halfX = p.first;
				if (p.second < cutoffY) {
					extr.push_back(pMax);
					break;
				}
			}
			p0.swap(p);
		}
	}
	if (!halfX || pMMax.second - pMin.second <= 4) {	// why 4? 5 maybe enough to identify a peak
		//_spec = eSpec::EVEN;
		Err(Spec(_spec) + SepSCl + sParams + " are not called").Throw(false);	// even distribution
	}
#ifdef MY_DEBUG
	cout << "pMMax: " << pMMax.first << TAB << pMMax.second << LF;
#endif
	//== define splined max point
	pMMax.second = 0;
	for (const auto& p : extr)
		if (p.second > pMMax.second)	pMMax = p;

	//== define if sequence is modulated
	bool isPeakAfterDip = false;
	{
		bool isDip = false;
		auto itv = extr.begin();		// always 0,0
		itv++;							// always 0,0 as well
		dVal_t val = (++itv)->second;	// first unzero value in sequence
		itv++;

		// from now odd is always dip, even - peak, last - peak
		// looking critical dip in extremes
		for (int i = 1; itv != extr.end(); val = itv++->second) {
			//if(i % 2)		cout << float(val - itv->second) / pMMax.second << "\tdip\n";
			//else			cout << float(itv->second - val) / pMMax.second << "\tpeak\n";
			if (i++ % 2)	// dip
				isDip = float(val - itv->second) / pMMax.second > 0.3;
			else 			// peak
				if (isPeakAfterDip = isDip && float(itv->second - val) / pMMax.second > 0.1)
					break;
		}
		if (isPeakAfterDip)	_spec = eSpec::MODUL;
	}

	//== set _base
	fraglen diffX = halfX - pMMax.first;
	_base = fraglen(float(diffX) * (isPeakAfterDip ? 0.9F : (diffX > 20 ? 0.1F : 0.35F)));
#ifdef MY_DEBUG
	cout << "isPeakAfterDip: " << isPeakAfterDip << "\thalfX: " << halfX << "\tdiffX: " << diffX << "\tbase: " << _base << LF;
#endif
}

void Distrib::PrintWarning(dostream& s)
{
	if (_base == FRAGLEN_MAX)	return;

	const string sInaccurate = " may be biased";

	if (_base == smoothBase)
		s << Spec(eSpec::SMOOTH) << LF;
	if (_summit.first - begin()->first<SSpliner<dVal_t>::SilentLength(eCurveType::SMOOTH, _base)
	|| begin()->second / _summit.second > 0.95)
		Err(Spec(eSpec::HTRIM) + SepSCl + sParams + sInaccurate).Warning();
	else if (_spec == eSpec::MODUL)
		s << Spec(_spec) << LF;
	else if (begin()->second / _summit.second > 0.5)
		Err(Spec(eSpec::TRIM)).Warning();
	else {
		ADPs::ADP adp;

		adp.CalcPCC(*this, _indADPs.GetBestIndex(), _summit.first, false);	// sorts params
		const float diffPCC = adp.PCC - _indADPs.GetBestPCC();

#ifdef MY_DEBUG
		s << "summit: " << summit.first << "\tPCCsummit: " << adp.PCC << "\tdiff PCC: " << diffPCC << LF;
#endif
		if (diffPCC > 0.01)
			Err(Spec(eSpec::DEFECT) + SepSCl + sParams + sInaccurate).Warning();
		else if (diffPCC > 0.002)
			Err(Spec(eSpec::SDEFECT)).Warning();
	}
}

void Distrib::PrintOriginal(dostream& dos) const
{
#ifdef _DUP_OUTPUT
	ofstream& s = dos.File();
#else
	auto& s = dos;
#endif
	const fraglen maxLen = INT_MAX / 10;

	s << "\nOriginal " << sDistrib << COLON << "\nlength\tfrequency\n";
	for (const value_type& f : *this) {
		if (f.first > maxLen)	break;
		s << f.first << TAB << f.second << LF;
	}
}


Distrib::Distrib(const char* fName, dostream& s)
{
	TabReader file(fName, FT::DIST);
	size_t cnt = 0;

	for (int x; file.GetNextLine();)
		if (x = file.UIntField(0))		// returns 0 if zero field is not an integer
			cnt += (*this)[x] = file.UIntField(1);
	if (cnt)
		s << SepCl << Size() << " records, " << cnt << " items";
}

void Distrib::ADParamsPrint(dostream& s, bool prWarning, bool prDistr)
{
	if (empty())
		s << "\nempty " << sDistrib << LF;
	else
		if (_base) {
			if (prWarning)	PrintWarning(s);
			_indADPs.Print(s);
			if (prDistr)	PrintOriginal(s);
		}
		else
			s << "\nDegenerate " << sDistrib << " (only " << size() << " points)\n";
	std::fflush(stdout);		// when called from a package
}

void Distrib::Print(dostream& s) const
{
	for (auto& f : *this)
		s << f.first << TAB << f.second << LF;
}