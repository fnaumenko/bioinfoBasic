/**********************************************************
Distrib.cpp
Last modified: 01/12/2025
***********************************************************/

#include "Distrib.h"
#include "spline.h"
#include <algorithm>    // std::sort
#include <utility>      // std::swap

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

static ADF ADFs[3] {
#define FACTOR1	cmFactors.first
#define FACTOR2	cmFactors.second
#define MEAN	p.first
#define SIGMA	p.second
#define SQR_SIGMA	SIGMA * SIGMA
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
	[](const fpair& p) ->fpair { return { SIGMA * SDPI, 2 * SQR_SIGMA}; },	// CommonFactors
	[](const fpair& p, fraglen x, const fpair& cmFactors) { 	// Value
		return exp(-pow((log(x) - MEAN), 2) / FACTOR2) / (FACTOR1 * x);
	},
	[](const fpair& p) { return exp(MEAN - SQR_SIGMA); },		// Mode
	[](const fpair& p) { return exp(MEAN + SQR_SIGMA / 2); },	// Mean
	[](const fpair& p) { return exp(MEAN); },					// Median
	[](const fpair& keypts, fpair& p) {							// CalcParams
		const float lgM = log(X_HIGH);		// logarifm of Mode
		const float lgH = log(X_HALFHIGH);	// logarifm of middle height
		MEAN = (lgM * (lghRatio + lgM - lgH) + (lgH * lgH - lgM * lgM) / 2) / lghRatio;
		SIGMA = sqrt(MEAN - lgM);
	},
	},
	{ "Gamma",
	[](const fpair& p) ->fpair { return { MEAN - 1, float(pow(SIGMA, MEAN)) }; },	// CommonFactors
	[](const fpair& p, fraglen x, const fpair& cmFactors) { // Value
		return pow(x, FACTOR1) * exp(-(x / SIGMA)) / FACTOR2;
	},
	[](const fpair& p) { return (MEAN - 1) * SIGMA; },		// Mode
	[](const fpair& p) { return MEAN * SIGMA; },			// Mean 
	[](const fpair& p) { return 0.f; },						// Median: undefined
	[](const fpair& keypts, fpair& p) {						// CalcParams
		SIGMA = (X_HALFHIGH - X_HIGH * (1 + log(X_HALFHIGH / X_HIGH))) / lghRatio;
		MEAN = (X_HIGH / SIGMA) + 1;
	}
	},
};

//===== Distrib::SetADParams

#define SETW left<<setw(4)
#define UNITAB SETW<<SPACE<<TAB	// tab stretching 4 spaces to display regardless of tab size (4 or 8)

void Distrib::SetADParams::IndADParams::Print(dostream& s, float maxPCC) const
{
	if (IsSet()) {
		auto& adfs = ADFs[Index];

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

bool Distrib::SetADParams::IsSetInSorted(eCType ctype) const
{
	const dind ind = GetDType(ctype);
	for (const auto& dp : _allParams)
		if (dp.Index == ind)
			return dp.IsSet();
	return false;
}

int Distrib::SetADParams::SetSortedCount() const
{
	int cnt = 0;
	for (const auto& dp : _allParams)
		cnt += dp.IsSet();
	return cnt;
}

void Distrib::SetADParams::Sort()
{
	if (!_sorted) {
		sort(_allParams.begin(), _allParams.end(),
			[](const ADParams& dp1, const ADParams& dp2) -> bool
			{ return dp1 > dp2; }
		);
		_sorted = true;
	}
}

Distrib::SetADParams::SetADParams()
{
	int i = 0;
	for (auto& dp : _allParams)
		dp.Index = i++;
}

void Distrib::SetADParams::Print(dostream& s)
{
	static const char* N[] = { "mean", "sigma" };	// normal, lognormal parameters
	static const char* G[] = { "alpha", "beta" };	// gamma parameters
	static const char* P[] = { "p1", "p2" };		// unified parameters
	static const char* a[] = { "* ", "**" };		// asterisks - footnotes
	const bool notSingle = SetSortedCount() > 1;	// more then 1 output distr type
	float maxPCC = 0;

	Sort();			// should already be sorted by PrintSpecs(), but just in case
	const bool isGamma = IsSetInSorted(eCType::GAMMA);

	// ** print title
	s << LF << UNITAB << " PCC\t";
	if (notSingle)
		s << "relPCC\t",
		maxPCC = _allParams[0].PCC;
	if (!isGamma)		s << N[0] << TAB << N[1];
	else if (notSingle)	s << P[0] << a[0] << TAB << P[1] << a[1];
	else				s << G[0] << TAB << G[1];
	s << "\tMode\tMean";
	if (notSingle || !isGamma)
		s << "\tMedian";

	// ** print values
	s << LF;
	for (const auto& params : _allParams)
		params.Print(s, maxPCC);

	// ** print note
	if (notSingle && isGamma) {
		auto title = ADFs[GetDType(eCType::GAMMA)].Title;
		s << LF;
		for (BYTE i = 0; i < 2; i++)
			s << setw(3) << a[i] << P[i] << " - " << N[i] << ", or "
			<< G[i] << " for " << title << LF;
	}
}


//===== Distrib

const char* Distrib::sDistrib = "distribution";

double Distrib::GetApprValue(eCType ctype, float mean, float sigma, fraglen x)
{
	const fpair p{ mean, sigma };
	return ADFs[GetDType(ctype)].GetValue(p, x);
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

fraglen Distrib::GetBase()
{
	using rpoint = pair<fraglen, dVal_t>;	// initial raw sequence point

	const int CutoffFrac = 100;	// fraction of the maximum height below which scanning stops on the first pass
	dVal_t	cutoffY = 0;		// Y-value below which scanning stops on the first pass
	fraglen halfX = 0;
	USHORT peakCnt = 0;
	bool up = false;
	auto it = begin();
	rpoint p0(*it), p;			// previous, current point
	rpoint pMin(0, 0), pMax(pMin), pMMax(pMin);	// current, previous, maximum point
	vector<value_type> extr;		// local extremums
	SSpliner<dVal_t> spliner(eCurveType::ROUGH, 1);

	//== define pMMax and halfX
	extr.reserve(20);
	for (it++; it != end(); p0 = p, it++) {
		p.first = it->first;
		p.second = dVal_t(spliner.Push(it->second));
		if (p.second > p0.second) {		// increasing part
			if (!up) {					// treat pit
				extr.push_back(pMax); pMin = p0; up = true;
			}
		}
		else {							// decreasing part
			if (up) {					// treat peak
				extr.push_back(pMin); pMax = p0; up = false;

				if (p0.second > pMMax.second) {
					pMMax = p0;
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
	}
	if (!halfX || pMMax.second - pMin.second <= 4) {	// why 4? 5 maybe enough to identify a peak
		//_spec = eSpec::EVEN;
		Err(Spec(_spec) + SepSCl + sParams + " are not called").Throw(false);	// even distribution
		return 0;
	}
#ifdef MY_DEBUG
	cout << "pMMax: " << pMMax.first << TAB << pMMax.second << LF;
#endif
	//== define splined max point
	pMMax.second = 0;
	for (auto& p : extr) {
		if (p.second > pMMax.second)	pMMax = p;
		//cout << p.first << TAB << p.second << LF;
	}

	//== define if sequence is modulated
	auto itv = extr.begin();	// always 0,0
	itv++;						// always 0,0 as well
	p0 = *(++itv);				// first point in sequence
	//pMin = pMMax;
	//pMax = make_pair(0, 0);
	int i = 1;
	bool isDip = false, isPeakAfterDip = false;

	// from now odd is always dip, even - peak, last - peak
	// looking critical dip in extr
	for (itv++; itv != extr.end(); i++, itv++) {
		p = *itv;	// we don't need point, using the Y-coordinate is enough. Point is used for debugging
		//cout << p.first << TAB << p.second << TAB;
		//if(i % 2)		cout << float(p0.second - p.second) / pMMax.second << "\tdip\n";
		//else			cout << float(p.second - p0.second) / pMMax.second << "\tpeak\n";
		if (i % 2)		// dip
			isDip = float(p0.second - p.second) / pMMax.second > 0.3;
		else 			// peak
			if (isPeakAfterDip = isDip && float(p.second - p0.second) / pMMax.second > 0.1)
				break;
		p0 = p;
	}
	// set  
	if (isPeakAfterDip)	_spec = eSpec::MODUL;

	if (!halfX)		return smoothBase;
	fraglen diffX = halfX - pMMax.first;
#ifdef MY_DEBUG
	fraglen base = fraglen((isPeakAfterDip ? 0.9F : (diffX > 20 ? 0.1F : 0.35F)) * diffX);
	cout << "isPeakAfterDip: " << isPeakAfterDip << "\thalfX: " << halfX << "\tdiffX: " << diffX << "\tbase: " << base << LF;
	return base;
#else
	return fraglen(float(diffX) * (isPeakAfterDip ? 0.9F : (diffX > 20 ? 0.1F : 0.35F)));
#endif
}

fpair Distrib::GetKeyPoints(fraglen base, dpoint& summit) const
{
	dpoint p0 = make_pair(begin()->first, float(begin()->second));	// previous point
	dpoint p{};					// current point
	SSpliner<dVal_t> spliner(
#ifdef MY_DEBUG					// to visualize SPIKED or SMOOTH distributions separately
		eCurveType::SPIKED, 
		//eCurveType::SMOOTH,
#else
		base <= smoothBase ? eCurveType::ROUGH : eCurveType::SMOOTH,
#endif
		base);

	summit.second = 0;
#ifdef MY_DEBUG
	_spline.clear();
	fpair keyPts(0, 0);
#endif
	for (const value_type& f : *this) {
		p.first = spliner.CorrectX(f.first);	// X: minus MA & MM base back shift
		p.second = spliner.Push(f.second);		// Y: splined
#ifdef MY_DEBUG
		if (_fillSpline)	_spline.push_back(p);	// to print
#endif
		if (p.second >= summit.second)	
			std::swap(p, summit);
		else {
			if (p.second < summit.second / hRatio) {
#ifdef MY_DEBUG
				if (!keyPts.first) {
					keyPts.first = float(summit.first);
					keyPts.second = p0.first + p0.second / (p.second + p0.second);
				}
			}
			if (p.second < summit.second / (hRatio * 5)) {
#endif
				break;
			}
			std::swap(p, p0);
		}
	}

#ifdef MY_DEBUG
	return keyPts;
#else
	return fpair(
		float(summit.first),							// summit X
		p0.first + p0.second / (p.second + p0.second)	// final point with half height; proportional X
	);
#endif
}

void Distrib::CalcPCC(dind ind, ADParams& dParams, fraglen Mode, bool full) const
{
	const fpair commFactors = ADFs[ind].CommonFactors(dParams.Params);	// two constant terms of the distrib equation
	const auto dVal = ADFs[ind].Value;	// function that calculates the 'type' distribution coordinate
	const double cutoffY = dVal(dParams.Params, Mode, commFactors) / 1000;	// break when Y became less then 0.1% of max value
	double	sumA = 0, sumA2 = 0;	// sum, sum of squares of original values
	double	sumB = 0, sumB2 = 0;	// sum, sum of squares of calculated values
	double	sumAB = 0;				// sum of products of original and calculated values
	UINT	cnt = 0;				// count of points

	// one pass PCC calculation
	dParams.SetUndefPcc();
	for (const value_type& f : *this) {
		if (!full && f.first < Mode)		continue;
		const double b = dVal(dParams.Params, f.first, commFactors);	// y-coordinate (value) of the calculated sequence
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
	if (!isNaN(pcc))	dParams.PCC = pcc;
}

void Distrib::CallParams(dind ind, fraglen base, dpoint& summit)
{
	auto calcParams = ADFs[ind].CalcParams;
	const BYTE failCntLim = 2;	// max count of base's decreasing steps after which PCC is considered only decreasing
	BYTE failCnt = 0;			// counter of base's decreasing steps after which PCC is considered only decreasing
	dpoint summit0;				// temporary summit
	ADParams dParams0, dParams;	// temporary, final  PCC & mean(alpha) & sigma(beta)
#ifdef MY_DEBUG
	int i = 0;					// counter of steps
#endif

	// calculate the highest PCC by iteratively searching through the 'base' values
	//base = 9;					// to print spline for fixed base
	//for (int i=0; !i; i++) {	// to print spline for fixed base
	for (; base; base--) {
		const fpair keypts = GetKeyPoints(base, summit0);

		calcParams(keypts, dParams0.Params);
		CalcPCC(ind, dParams0, summit0.first);
#ifdef MY_DEBUG
		* _s << setw(4) << setfill(SPACE) << left << ++i;
		*_s << "base: " << setw(2) << base << "  summitX: " << keypts.first << "\tpcc: " << dParams0.PCC;
		if (dParams0 > dParams)	*_s << "\t>";
		*_s << LF;
		if (_fillSpline) { for (dpoint p : _spline)	*_s << p.first << TAB << p.second << LF; _fillSpline = false; }
#endif
		if (dParams0 > dParams) {
			dParams = dParams0;
			summit = summit0;
			failCnt = 0;
		}
		else {
			if (dParams0.PCC > 0)	failCnt++;		// negative PCC is possible in rare cases
			else if (dParams0.IsUndefPcc()) {
				dParams.SetUndefPcc();
				break;
			}
			if (failCnt > failCntLim)	break;
		}
	}
	_allParams.SetParams(ind, dParams);

#ifdef MY_DEBUG
	* _s << LF;
#endif
}

void Distrib::PrintSpecs(dostream& s, fraglen base, const Distrib::dpoint& summit)
{
	const string sInaccurate = " may be biased";

	if (base == smoothBase)
		s << Spec(eSpec::SMOOTH) << LF;
	if (summit.first - begin()->first<SSpliner<dVal_t>::SilentLength(eCurveType::SMOOTH, base)
	|| begin()->second / summit.second > 0.95)
		Err(Spec(eSpec::HTRIM) + SepSCl + sParams + sInaccurate).Warning();
	else if (_spec == eSpec::MODUL)
		s << Spec(_spec) << LF;
	else if (begin()->second / summit.second > 0.5)
		Err(Spec(eSpec::TRIM)).Warning();
	else {
		ADParams dParams;

		CalcPCC(_allParams.GetBestIndex(), dParams, summit.first, false);	// sorts params
		const float diffPCC = dParams.PCC - _allParams.GetBestPCC();

#ifdef MY_DEBUG
		s << "summit: " << summit.first << "\tPCCsummit: " << dParams.GetBestPCC() << "\tdiff PCC: " << diffPCC << LF;
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

//#define _TIME
#ifdef _TIME
#include <chrono> 
using namespace std::chrono;
#endif

void Distrib::Print(dostream& s, eCType ctype, bool prWarning, bool prDistr)
{
	if (empty())		s << "\nempty " << sDistrib << LF;
	else {
		fraglen base = GetBase();	// initialized returned value
		if (base) {
#ifdef _TIME
			auto start = high_resolution_clock::now();
			const int	tmCycleCnt = 1000;
#endif			
			// For optimization purposes, we can initialize base, keypts & summit at the first call of CallParams,
			// and use them on subsequent calls to avoid repeated PCC iterations.
			// However, the same base (and, as a consequence, keypts & summit) only works well for LNORM and GAMMA.
			// For the best NORM, base may be less, therefore, for simplicity and reliability, all parameters are always recalculated
#ifdef MY_DEBUG
			_s = &s;
			if (_fillSpline)	_spline.reserve(size() / 2);
#endif
#ifdef _TIME
			for (int i = 0; i < tmCycleCnt; i++)
#endif
				dpoint summit;				// returned value
			for (dind i = 0; i < eCType::CNT; i++)
				if (IsIndex(ctype, i))
					CallParams(i, base, summit);
#ifdef _TIME
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			s << duration.count() / tmCycleCnt << " mcs\n";
#else
			// check for NORM if LNORM is defined
			if (IsType(ctype, eCType::LNORM) && !IsType(ctype, eCType::NORM)) {
				CallParams(GetDType(eCType::NORM), base, summit);
				_allParams.ClearNormDistBelowThreshold(1.02F);	// threshold 2%
			}
#endif
			if (prWarning)	PrintSpecs(s, base, summit);
			_allParams.Print(s);
			if (prDistr)	PrintOriginal(s);
		}
		else
			s << "\nDegenerate " << sDistrib << " (only " << size() << " points)\n";
	}
	std::fflush(stdout);		// when called from a package
}
