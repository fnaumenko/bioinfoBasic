/**********************************************************
Distrib.h
2023 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 01/14/2025
-------------------------
Provides value (typically frequency) distribution functionality
***********************************************************/
#pragma once

#include "DataReader.h"
#include <array>

// MY_DEBUG should be managed via DataReader.h

using dVal_t = size_t;	// type of distribution value

// 'Distrib' represents a value frequency distribution and its approximation by a two-parameter distribution
class Distrib : map<int, dVal_t>	// sync type with 'rpoint' in .cpp 
{
public:
	// combined type of distribution
	enum /*class*/ eDType {		// not class to have a cast to integer by default
		NORM = 1 << 0,
		LNORM = 1 << 1,
		GAMMA = 1 << 2,
		CNT = 3,
	};

	// method of smoothing distribution
	//enum eSmooth {
	//	SPLINE,		// sliding splining
	//	INTERPOL,	// Bezier interpolation
	//};

	static const char* sDistrib;

	// Default constructor
	Distrib() {}

	// Constructor by ready distribution file
	//	@param fname: name of ready distribution file
	Distrib(const char* fname, dostream& s);

	// Returns size distribution
	size_t Size() const { return size(); }

	// Increments value frequency
	void IncrFreq(fraglen val) { (*this)[val]++; }

	// Returns the value the value of the approximate distribution function at a given point
	//	@param dtype: type of distribution
	//	@param mean: mean (for norm, lognorm) or alpha (for gamma)
	//	@param sigma: sigma (for norm, lognorm) or beta (for gamma)
	//	@param x: X-value of the point
	//	@returs Y-value of the point
	static double GetApprValue(eDType dtype, float mean, float sigma, fraglen x);

	// Calculate approximate distribution parameters
	//	@param type: combined type of distribution
	//	@param smooth: method of smoothing
	void CalcADParams(eDType type/*, eSmooth smooth = eSmooth::SPLINE*/);

	// Prints approximate distribution parameters on a new line
	//	@param s[out]: print stream
	//	@param prWarning[in]: if true then print possible warning message
	//	@param prDistr[in]: if true then print original distribution additionally
	void Print(dostream& s, bool prWarning, bool prDistr);

private:
	using dind = BYTE;						// inner distribution index
	using dpoint = pair<fraglen, float>;	// distribution point 

	// Returns combined distribution type by inner distribution index
	static eDType GetCType(dind ind) { return eDType(1 << ind); }

	// Returns inner distribution index by combined distribution type
	const static dind GetDType(eDType dtype) { return RightOnePos(int(dtype)); }

	enum class eSpec {	// distribution specification
		CLEAR,		// normal quality;	exclusive
		SMOOTH,		// complementary
		MODUL,		// modulated;	complementary
		EVEN,		// flat;	exclusive
		TRIM,		// trimmed;	exclusive
		HTRIM,		// heavily trimmed;		exclusive
		SDEFECT,	// slightly defective;	exclusive
		DEFECT		// defective; exclusive
	};

	const fraglen smoothBase = 1;	// splining base for the smooth distribution

	// Keeps approximate distribution parameters: PCC, mean(alpha), sigma(beta)
	struct ADParams
	{
		float	PCC = 0;	// Pearson correlation coefficient
		fpair	Params{};	// mean(alpha), sigma(beta)

		bool operator >(const ADParams& dp) const { return PCC > dp.PCC; }

		bool IsUndefPcc() const { return PCC == -1; };

		void SetUndefPcc() { PCC = -1; };
	};

	// 'SetADParams' represents a collection of approximate distribution parameters for all type of distribution
	class SetADParams
	{
		// Indexed ADParams: struct ADParams supplied with inner index
		struct IndADParams : public ADParams
		{
			dind	Index;		// inner index

			// Returns true if AD parameters set
			bool IsSet() const { return PCC; }

			void Copy(const ADParams& dp) { PCC = dp.PCC; Params = dp.Params; }

			// Prints AD parameters
			//	@param s: print stream
			//	@param maxPCC: masimum PCC to print relative PCC percentage
			void Print(dostream& s, float maxPCC) const;
		};

		array<IndADParams, eDType::CNT>	_setADParams;
		bool _sorted = false;

		// Returns true if AD parameters set in sorted instance
		bool IsSetInSorted(eDType dtype) const;

		// Returns number of AD parameters set in sorted instance
		int SetSortedCount() const;

		// Returns AD Params by combined distribution type
		ADParams& Params(eDType dtype) { return _setADParams[GetDType(dtype)]; }

		// Sorts in PCC descending order
		void Sort();

	public:
		// Default constructor
		SetADParams();

		float GetBestPCC() const { return _setADParams[0].PCC; }

		// Calculates and set approximate distribution parameters by index
		//	@param ind: inner distribution index
		//	@param adp: approximate distribution parameters
		void SetParams(dind ind, const ADParams& adp) { _setADParams[ind].Copy(adp); }

		// Clear normal distribution if its PCC is less then lognorm PCC by the threshold
		void ClearNormDistBelowThreshold(float thresh) {
			if (Params(eDType::LNORM).PCC / Params(eDType::NORM).PCC > thresh)
				Params(eDType::NORM).PCC = 0;
		}

		// Sorts parameters and returns inner index of distribution with the highest PCC
		//	@returns inner index of distribution with the highest (best) PCC
		dind GetBestIndex() { Sort(); return _setADParams[0].Index; }

		// Prints sorted distibutions params on a new line
		//	@param s: output stream
		void Print(dostream& s);
	};

	// Returns specification string by specification type
	static const string Spec(eSpec s);

	// Returns true if inner index is represented in combo cType
	static bool IsIndex(eDType cType, dind ind) { return cType & (1 << ind); }

	// Returns true if exclusive type is represented in combo cType
	//	@param test: test combo cType
	//	@param excl: exclusive cType
	static bool IsType(eDType test, eDType excl) { return test & excl; }

	eSpec _spec = eSpec::CLEAR;		// distribution specification
	SetADParams	_setADParams;		// approximate distribution parameters for all type of distributions
	// these two fields are needed to print warnings after calculating the parameters
	fraglen	_base = FRAGLEN_MAX;	// moving window half-length of best spline
	dpoint	_summit;				// X,Y coordinates of best splined (smoothed) summit

	//eSmooth	_smooth = eSmooth::SPLINE;
#ifdef MY_DEBUG
	mutable vector<dpoint> _spline;		// splining curve (container) to visualize splining
	mutable bool _fillSpline = true;	// true if fill splining curve (container)
	dostream* _s = NULL;				// print stream
#endif

	// Set moving window half-length of appropriate spline (estimated base)
	void SetBase();

	// Builds spline curve and defines key points
	//	@param base[in]: moving window half-length
	//	@param summit[out]: returned X,Y coordinates of splined (smoothed) summit
	//	@returns key points: X-coord of highest point, X-coord of right middle hight point
	fpair GetKeyPoints(fraglen base, dpoint& summit) const;

	// Compares this sequence with calculated one with given mean&sigma, and returns PCC
	//	@param ind[in]: inner distribution index
	//	@param dParams[in, out]: returned PCC, input mean(alpha) & sigma(beta)
	//	@param Mode[in]: X-coordinate of summit
	//	@param full[in]: if true then correlate from the beginning, otherwiase from summit
	//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
	void CalcPCC(dind ind, ADParams& dParams, int Mode, bool full = true) const;

	// Calculates the best approximate distribution parameters for a specific type of distribution
	//	@param ind[in]: inner distribution index
	void SetParamsForSpline(dind ind);

	// Calculates the best approximate distribution parameters for a specific type of distribution
	//	@param ind[in]: inner distribution index
	//void SetParamsForInterpol(dind ind) {}

	// Calculates the best approximate distribution parameters for a specific type of distribution
	//	@param ind[in]: inner distribution index
	void SetParams(dind ind) { SetParamsForSpline(ind); }

	// Prints warnings about poor quality distributions
	//	@param s: print stream
	void PrintWarning(dostream& s);

	// Prints original distribution as a set of <value>-<size> pairs
	//	@param s: print stream
	void PrintOriginal(dostream& s) const;
};
