/**********************************************************
Distrib.h
2023 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 01/11/2025
-------------------------
Provides value (typically frequency) distribution functionality
***********************************************************/
#pragma once

#include "DataReader.h"
#include <array>

//#define MY_DEBUG

using dVal_t = size_t;	// type of distribution value

// 'Distrib' represents a value frequency distribution and its approximation by a two-parameter distribution
class Distrib : map<fraglen, dVal_t>
{
public:
	// combined type of distribution
	enum /*class*/ eCType {		// not class to have a cast to integer by default
		NORM = 1 << 0,
		LNORM = 1 << 1,
		GAMMA = 1 << 2,
		CNT = 3,
	};
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
	//	@param ctype: type of distribution
	//	@param mean: mean (for norm, lognorm) or alpha (for gamma)
	//	@param sigma: sigma (for norm, lognorm) or beta (for gamma)
	//	@param x: X-value of the point
	//	@returs Y-value of the point
	static double GetApprValue(eCType ctype, float mean, float sigma, fraglen x);

	// Calculate and print distribution on a new line
	//	@param s[out]: print stream
	//	@param type[in]: combined type of distribution
	//	@param prWarning[in]: if true then print possible warning message
	//	@param prDistr[in]: if true then print original distribution additionally
	void Print(dostream& s, eCType type, bool prWarning, bool prDistr);

private:
	using dind = BYTE;						// inner distribution index
	using dpoint = pair<fraglen, float>;	// distribution point 

	// Returns combined distribution type by inner distribution index
	static eCType GetCType(dind ind) { return eCType(1 << ind); }

	// Returns inner distribution index by combined distribution type
	const static dind GetDType(eCType ctype) { return RightOnePos(int(ctype)); }

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

	static const char* sTitle[];
	static const string sSpec[];
	static const string sParams;
	static const string sInaccurate;
	const fraglen smoothBase = 1;	// splining base for the smooth distribution

	// Keeps approximation distribution parameters: PCC, mean(alpha), sigma(beta)
	struct ADParams
	{
	private:
		static const float UndefPCC;
	public:
		float	PCC = 0;		// Pearson correlation coefficient
		fpair	Params{};		// mean(alpha), sigma(beta)

		bool operator >(const ADParams& dp) const { return PCC > dp.PCC; }

		bool IsUndefPcc() const { return PCC == UndefPCC; };

		void SetUndefPcc() { PCC = UndefPCC; };
	};

	// 'SetADParams' represents a collection of approximation distribution parameters for all type of distribution
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

		array<IndADParams, eCType::CNT>	_allParams;
		bool _sorted = false;

		// Returns true if AD parameters set in sorted instance
		bool IsSetInSorted(eCType ctype) const;

		// Returns number of AD parameters set in sorted instance
		int SetSortedCount() const;

		// Returns AD Params by combined distribution type
		ADParams& Params(eCType ctype) { return _allParams[GetDType(ctype)]; }

		// Sorts in PCC descending order
		void Sort();

	public:
		// Default constructor
		SetADParams();

		float GetBestPCC() const { return _allParams[0].PCC; }

		// Set distribution parameters by index
		//	@param ind: inner distribution index
		//	@param adp: approximation distribution parameters
		void SetParams(dind ind, const ADParams& adp) { _allParams[ind].Copy(adp); }

		// Clear normal distribution if its PCC is less then lognorm PCC by the threshold
		void ClearNormDistBelowThreshold(float thresh) {
			if (Params(eCType::LNORM).PCC / Params(eCType::NORM).PCC > thresh)
				Params(eCType::NORM).PCC = 0;
		}

		// Sorts parameters and returns inner index of distribution with the highest PCC
		//	@returns inner index of distribution with the highest (best) PCC
		dind GetBestIndex() { Sort(); return _allParams[0].Index; }

		// Prints sorted distibutions params on a new line
		//	@param s: output stream
		void Print(dostream& s);
	};

	// Returns specification string by specification type
	static const string Spec(eSpec s) { return "Distribution " + sSpec[int(s)]; }

	// Returns true if inner index is represented in combo cType
	static bool IsIndex(eCType cType, dind ind) { return cType & (1 << ind); }

	// Returns true if combo type is represented in combo cType
	static bool IsType(eCType cType, eCType type) { return cType & type; }

	eSpec _spec = eSpec::CLEAR;		// distribution specification
	SetADParams	_allParams;			// distributions parameters
#ifdef MY_DEBUG
	mutable vector<dpoint> _spline;		// splining curve (container) to visualize splining
	mutable bool _fillSpline = true;	// true if fill splining curve (container)
	dostream* _s = NULL;				// print stream
#endif

	// Returns estimated moving window half-length ("base")
	//	@returns estimated base, or 0 in case of degenerate distribution
	fraglen GetBase();

	// Builds spline curve and defines key points
	//	@param base: moving window half-length
	//	@param summit: returned X,Y coordinates of spliced (smoothed) summit
	//	@returns key points: X-coord of highest point, X-coord of right middle hight point
	fpair GetKeyPoints(fraglen base, dpoint& summit) const;

	// Compares this sequence with calculated one with given mean&sigma, and returns PCC
	//	@param ind[in]: inner distribution index
	//	@param dParams[in, out]: returned PCC, input mean(alpha) & sigma(beta)
	//	@param Mode[in]: X-coordinate of summit
	//	@param full[in]: if true then correlate from the beginning, otherwiase from summit
	//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
	void CalcPCC(dind ind, ADParams& dParams, fraglen Mode, bool full = true) const;

	// Calculates the best distribution parameters
	//	@param ind[in]: inner distribution index
	//	@param base[in]: moving window half-length
	//	@param summit[out]: returned X,Y coordinates of best spliced (smoothed) summit
	void CallParams(dind ind, fraglen base, dpoint& summit);

	// Prints original distribution specification (flaws)
	//	@param s: print stream
	//	@param base: moving window half-length
	//	@param summit: X,Y coordinates of spliced (smoothed) summit
	void PrintSpecs(dostream& s, fraglen base, const dpoint& summit);

	// Prints original distribution as a set of <value>-<size> pairs
	//	@param s: print stream
	void PrintOriginal(dostream& s) const;
};
