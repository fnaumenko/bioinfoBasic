/**********************************************************
Distrib.cpp
Last modified: 04/29/2024
***********************************************************/

#include "Distrib.h"
#include "spline.h"
#include <algorithm>    // std::sort

const float SDPI = float(sqrt(3.1415926 * 2));		// square of doubled Pi
// ratio of the summit height to height of the measuring point
const int hRatio = 2;
// log of ratio of the summit height to height of the measuring point
const float lghRatio = float(log(hRatio));

const char* Distrib::sDistrib = "distribution";
const char* Distrib::sTitle[] = { "Norm", "Lognorm", "Gamma" };
const float Distrib::DParams::UndefPCC = -1;
const string Distrib::sParams = "parameters";
const string Distrib::sInaccurate = " may be inaccurate";
const string Distrib::sSpec[] = {
	"is degenerate",
	"is smooth",
	"is modulated",
	"is even",
	"is cropped to the left",
	"is heavily cropped to the left",
	"looks slightly defective on the left",
	"looks defective on the left"
};

// Returns two constant terms of the distrib equation of type, supplied as an index
//	@param p: distrib params: mean/alpha and sigma/beta
//	@returns: two initialized constant terms of the distrib equation
//	Implemeted as an array of lambdas treated like a regular function and assigned to a function pointer.
fpair(*GetEqTerms[])(const fpair& p) = {
	[](const fpair& p) -> fpair { return { p.second * SDPI, 0.f}; },						// normal
	[](const fpair& p) -> fpair { return { p.second * SDPI, 2 * p.second * p.second}; },	// lognormal
	[](const fpair& p) -> fpair { return { p.first - 1, float(pow(p.second, p.first)) }; }	// gamma
};

// Returns y-coordinate by x-coordinate of the distrib of type, supplied as an index
//	@param p: distrib params: mean/alpha and sigma/beta
//	@param x: x-coordinate
//	@param eqTerms: two constant terms of the distrib equation
double (*Distrs[])(const fpair& p, fraglen x, const fpair& eqTerms) = {
	[](const fpair& p, fraglen x, const fpair& eqTerms) ->	double { return 		// normal
		exp(-pow(((x - p.first) / p.second), 2) / 2) / eqTerms.first; },
	[](const fpair& p, fraglen x, const fpair& eqTerms) ->	double { return 		// lognormal
		exp(-pow((log(x) - p.first), 2) / eqTerms.second) / (eqTerms.first * x); },
	[](const fpair& p, fraglen x, const fpair& eqTerms) ->	double { return 		// gamma
		pow(x, eqTerms.first) * exp(-(x / p.second)) / eqTerms.second; }
};

double Distrib::GetVal(eCType ctype, float mean, float sigma, fraglen x)
{
	const fpair p{ mean, sigma };
	auto type = GetDType(ctype);

	return Distrs[type](p, x, GetEqTerms[type](p));
}

// Returns distribution mode
//	@param p: distrib params: mean/alpha and sigma/beta
float (*GetMode[])(const fpair& p) = {
	[](const fpair& p) { return 0.0f; },								// normal
	[](const fpair& p) { return exp(p.first - p.second * p.second); },	// lognormal
	[](const fpair& p) { return (p.first - 1) * p.second; }				// gamma 
};

// Returns distribution Mean (expected value)
//	@param p: distrib params: mean/alpha and sigma/beta
float (*GetMean[])(const fpair& p) = {
	[](const fpair& p) { return 0.0f; },									// normal
	[](const fpair& p) { return exp(p.first + p.second * p.second / 2); },	// lognormal
	[](const fpair& p) { return p.first * p.second; }						// gamma 
};

// Calculates distribution parameters
//	@param keypts[in]: key points: X-coord of highest point, X-coord of right middle hight point
//	@param p[out]: returned params: mean(alpha) & sigma(beta)
void (*CalcParams[])(const fpair& keypts, fpair& p) = {
	[](const fpair& keypts, fpair& p) {		//*** normal
		p.first = keypts.first;
		p.second = float(sqrt(pow(keypts.second - p.first, 2) / lghRatio / 2));
	},
	[](const fpair& keypts, fpair& p) {		//*** lognormal
		const float lgM = log(keypts.first);		// logarifm of Mode
		const float lgH = log(keypts.second);		// logarifm of middle height

		p.first = (lgM * (lghRatio + lgM - lgH) + (lgH * lgH - lgM * lgM) / 2) / lghRatio;
		p.second = sqrt(p.first - lgM);
	},
	[](const fpair& keypts, fpair& p) {		//*** gamma
		p.second = (keypts.second - keypts.first * (1 + log(keypts.second / keypts.first))) / lghRatio;
		p.first = (keypts.first / p.second) + 1;
	}
};

void Distrib::AllDParams::QualDParams::Print(dostream& s, float maxPCC) const
{
	if (IsSet()) {
		s << Title << TAB;
		if (dParams.IsUndefPcc())
			s << "parameters cannot be called";
		else {
			s << setprecision(5) << dParams.PCC << TAB;
			if (maxPCC) {		// print percent to max PCC
				if (maxPCC != dParams.PCC)	s << setprecision(3) << 100 * ((dParams.PCC - maxPCC) / maxPCC) << '%';
				s << TAB;
			}
			s << setprecision(4) << dParams.Params.first << TAB << dParams.Params.second << TAB;
			const dtype type = GetDType(Type);
			float Mode = GetMode[type](dParams.Params);
			if (Mode)		// for normal Mode is equal to 0
				s << Mode << TAB << GetMean[type](dParams.Params);
		}
		s << LF;
	}
}

bool Distrib::AllDParams::IsSetInSorted(eCType ctype) const
{
	for (const QualDParams& dp : _allParams)
		if (dp.Type == ctype)
			return dp.IsSet();
	return false;
}

int Distrib::AllDParams::SetCntInSorted() const
{
	int cnt = 0;
	for (const QualDParams& dp : _allParams)
		cnt += dp.IsSet();
	return cnt;
}

void Distrib::AllDParams::Sort()
{
	if (!_sorted) {
		sort(_allParams.begin(), _allParams.end(),
			[](const QualDParams& dp1, const QualDParams& dp2) -> bool
			{ return dp1.dParams > dp2.dParams; }
		);
		_sorted = true;
	}
}

Distrib::AllDParams::AllDParams()
{
	int i = 0;
	for (QualDParams& dp : _allParams)
		dp.SetTitle(i++);
	//for (int i = 0; i < DTCNT; i++)
	//	_allParams[i].SetTitle(i);
}

Distrib::dtype Distrib::AllDParams::GetBestParams(DParams& dParams)
{
	Sort();
	const QualDParams& QualDParams = _allParams[0];
	dParams = QualDParams.dParams;
	return GetDType(QualDParams.Type);
}

void Distrib::AllDParams::Print(dostream& s)
{
	static const char* N[] = { "mean", "sigma" };
	static const char* G[] = { "alpha", "beta" };
	static const char* P[] = { "p1", "p2" };
	static const char* a[] = { "*", "**" };
	const bool notSingle = SetCntInSorted() > 1;	// more then 1 output distr type
	float maxPCC = 0;
	bool note = false;

	Sort();			// should already be sorted by PrintTraits(), but just in case
	s << LF << "\t PCC\t";
	if (notSingle)
		s << "relPCC\t",
		maxPCC = _allParams[0].dParams.PCC;
	if (!IsSetInSorted(eCType::GAMMA))	s << N[0] << TAB << N[1];
	else if (note = notSingle)			s << P[0] << a[0] << TAB << P[1] << a[1];
	else								s << G[0] << TAB << G[1];
	if (notSingle || !IsSetInSorted(eCType::NORM))
		s << "\tmode\texp.val";
	s << LF;
	for (const QualDParams& params : _allParams)
		params.Print(s, maxPCC);
	if (note) {
		s << LF;
		for (int i = 0; i < 2; i++)
			s << setw(3) << a[i] << P[i] << " - " << N[i] << ", or "
			<< G[i] << " for " << Distrib::sTitle[GetDType(eCType::GAMMA)] << LF;
	}
}

fraglen Distrib::GetBase()
{
	const int CutoffFrac = 100;	// fraction of the maximum height below which scanning stops on the first pass
	dVal_t	cutoffY = 0;		// Y-value below which scanning stops on the first pass
	fraglen halfX = 0;
	USHORT peakCnt = 0;
	bool up = false;
	auto it = begin();
	spoint p0(*it), p;			// previous, current point
	spoint pMin(0, 0), pMax(pMin), pMMax(pMin);	// current, previous, maximum point
	vector<spoint> extr;		// local extremums
	SSpliner<dVal_t> spliner(eCurveType::SPIKED, 1);

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
	pMMax = make_pair(0, 0);
	for (spoint p : extr) {
		if (p.second > pMMax.second)	pMMax = p;
		//cout << p.first << TAB << p.second << LF;
	}

	//== define if sequence is modulated
	auto itv = extr.begin();	// always 0,0
	itv++;						// always 0,0 as well
	p0 = *(++itv);				// first point in sequence
	pMin = pMMax;
	pMax = make_pair(0, 0);
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
#ifdef MY_DEBUG					// to visualize SPIKED or SMOOTH distributions individually
		eCurveType::SPIKED, 
		//eCurveType::SMOOTH,
#else
		base <= smoothBase ? eCurveType::SPIKED : eCurveType::SMOOTH,
#endif
		base);

	summit.second = 0;
#ifdef MY_DEBUG
	_spline.clear();
	fpair keyPts(0, 0);
#endif
	for (const value_type& f : *this) {
		p.first = spliner.CorrectX(f.first);		// X: minus MA & MM base back shift
		p.second = spliner.Push(f.second);	// Y: splined
#ifdef MY_DEBUG
		if (_fillSpline)	_spline.push_back(p);	// to print
#endif
		if (p.second >= summit.second)	summit = p;
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
			p0 = p;
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

void Distrib::CalcPCC(dtype type, DParams& dParams, fraglen Mode, bool full) const
{
	const auto dist = Distrs[type];								// function that calculates the 'type' distribution coordinate
	const fpair eqTerms = GetEqTerms[type](dParams.Params);	// two constant terms of the distrib equation
	const double cutoffY = dist(dParams.Params, Mode, eqTerms) / 1000;	// break when Y became less then 0.1% of max value
	double	sumA = 0, sumA2 = 0;	// sum, sum of squares of original values
	double	sumB = 0, sumB2 = 0;	// sum, sum of squares of calculated values
	double	sumAB = 0;				// sum of products of original and calculated values
	UINT	cnt = 0;				// count of points

	// one pass PCC calculation
	dParams.SetUndefPcc();
	for (const value_type& f : *this) {
		if (!full && f.first < Mode)		continue;
		const double b = dist(dParams.Params, f.first, eqTerms);	// y-coordinate (value) of the calculated sequence
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

void Distrib::CallParams(dtype type, fraglen base, dpoint& summit)
{
	const BYTE failCntLim = 2;	// max count of base's decreasing steps after which PCC is considered only decreasing
	BYTE failCnt = 0;			// counter of base's decreasing steps after which PCC is considered only decreasing
	dpoint summit0;				// temporary summit
	DParams dParams0, dParams;	// temporary, final  PCC & mean(alpha) & sigma(beta)
#ifdef MY_DEBUG
	int i = 0;					// counter of steps
#endif

	// calculate the highest PCC by iteratively searching through the 'base' values
	//base = 9;					// to print spline for fixed base
	//for (int i=0; !i; i++) {	// to print spline for fixed base
	for (; base; base--) {
		const fpair keypts = GetKeyPoints(base, summit0);

		CalcParams[type](keypts, dParams0.Params);
		CalcPCC(type, dParams0, summit0.first);
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
	_allParams.SetParams(type, dParams);

#ifdef MY_DEBUG
	* _s << LF;
#endif
}

void Distrib::PrintTraits(dostream& s, fraglen base, const Distrib::dpoint& summit)
{
	if (base == smoothBase)		s << Spec(eSpec::SMOOTH) << LF;
	if (summit.first - begin()->first < SSpliner<dVal_t>::SilentLength(eCurveType::SMOOTH, base)
	|| begin()->second / summit.second > 0.95)
		Err(Spec(eSpec::HCROP) + SepSCl + sParams + sInaccurate).Warning();
	else if (_spec == eSpec::MODUL)
		s << Spec(_spec) << LF;
	else if (begin()->second / summit.second > 0.5)
		Err(Spec(eSpec::CROP)).Warning();
	else {
		DParams dParams, dParams0;
		CalcPCC(_allParams.GetBestParams(dParams), dParams0, summit.first, false);	// sorts params
		float diffPCC = dParams0.PCC - dParams.PCC;
#ifdef MY_DEBUG
		s << "summit: " << summit.first << "\tPCCsummit: " << dParams.PCC << "\tdiff PCC: " << diffPCC << LF;
#endif
		if (diffPCC > 0.01)
			Err(Spec(eSpec::DEFECT) + SepSCl + sParams + sInaccurate).Warning();
		else if (diffPCC > 0.002)
			Err(Spec(eSpec::SDEFECT)).Warning();
	}
}

void Distrib::PrintSeq(ofstream& s) const
{
	const fraglen maxLen = INT_MAX / 10;

	s << "\nOriginal " << sDistrib << COLON << "\nlength\tfrequency\n";
	for (const value_type& f : *this) {
		if (f.first > maxLen)	break;
		s << f.first << TAB << f.second << LF;
	}
}

Distrib::Distrib(const char* fName)
{
	TabReader file(fName, FT::DIST);

	for (int x; file.GetNextLine();)
		if (x = file.UIntField(0))	// returns 0 if zero field is not an integer
			(*this)[x] = file.UIntField(1);
}

//#define _TIME
#ifdef _TIME
#include <chrono> 
using namespace std::chrono;
#endif

void Distrib::Print(dostream& s, eCType ctype, bool prDistr)
{
	if (empty())		s << "\nempty " << sDistrib << LF;
	else {
		fraglen base = GetBase();	// initialized returned value
		//if (IsDegenerate())
		if (!base)
			s << "\nDegenerate " << sDistrib << " (only " << size() << " points)\n";
		else {
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
			for (dtype i = 0; i < eCType::CNT; i++)
				if (IsType(ctype, i))
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
			PrintTraits(s, base, summit);
			_allParams.Print(s);
		}
		if (prDistr)	PrintSeq(s.File());
	}
	fflush(stdout);		// when called from a package
}
