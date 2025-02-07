/**********************************************************
Distrib.h
2023 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 02/07/2025
-------------------------
Provides value (typically frequency) distribution functionality
***********************************************************/
#pragma once

#include "common.h"
#include <map>
#include <array>

// MY_DEBUG should be managed via DataReader.h

using dVal_t = chrlen;	// type of distribution value
using dmap = map<int, dVal_t>;	// distribution map

// 'Distrib' represents a value frequency distribution and its approximation by a two-parameter distribution
class Distrib : dmap	// sync type with 'rpoint' in  .cpp 
{
public:
	// combined type of distribution
	enum eDType {		// not class to have a cast to integer by default
		NORM = 1 << 0,
		LNORM = 1 << 1,
		GAMMA = 1 << 2,
		CNT = 3,
	};

	// distribution smoothing method
	enum eDSmooth {
		SPLINE,		// sliding splining
		INTERPOL,	// Bezier interpolation
	};

	static const char* sDistrib;
	static const fraglen smoothBase = 1;	// splining base for the smooth distribution

	// Default constructor
	Distrib() {}

	// Constructor by ready distribution file
	//	@param fname: name of ready distribution file
	Distrib(const char* fname, dostream& s);

	// Returns size distribution
	size_t Size() const { return size(); }

	// Increments value frequency
	void IncrFreq(int val) { (*this)[val]++; }

	// Returns the value the value of the approximate distribution function at a given point
	//	@param dtype: type of distribution
	//	@param mean: mean (for norm, lognorm) or alpha (for gamma)
	//	@param sigma: sigma (for norm, lognorm) or beta (for gamma)
	//	@param x: X-value of the point
	//	@returs Y-value of the point
	static double GetApprValue(eDType dtype, float mean, float sigma, fraglen x);

	fpair	GetADParams() { return _indADPs.GetBestParams(); }

	// Calculate approximate distribution parameters
	//	@param type: combined type of distribution
	//	@param smooth: method of smoothing
	void CalcADParams(eDType type, eDSmooth smooth = eDSmooth::SPLINE) {
		if (!empty())	_indADPs.CalcParams(type, smooth, *this);
	}

	// Prints approximate distribution parameters on a new line
	//	@param s[out]: print stream
	//	@param prWarning[in]: if true then print possible warning message
	//	@param prDistr[in]: if true then print original distribution additionally
	void ADParamsPrint(dostream& s, bool prWarning, bool prDistr);

	// Prints original sequence as a set of lines of format "X-coord<TAB>Y-value"
	//	@param s[out]: print stream
	void Print(dostream& s) const;

private:
	using dpoint = pair<fraglen, float>;	// distribution point 

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

	// 'ADPs' represents a collection of approximate distribution parameters for all type of distribution
	class ADPs
	{
	public:
		// indexed type of distribution
		enum eDIndex {
			iNORM,
			iLNORM,
			iGAMMA,
		};

		// 'ADP' keeps Approximate Distribution Parameters: PCC, mean(alpha), sigma(beta)
		// The class is nested in ADPs only because of the use of the common 'eDIndex' type
		struct ADP
		{
			float	PCC = 0;	// Pearson correlation coefficient
			fpair	Params{};	// mean(alpha), sigma(beta)

			bool operator >(const ADP& dp) const { return PCC > dp.PCC; }

			bool IsUndefPcc() const { return PCC == -1; };

			void SetUndefPcc() { PCC = -1; };

			// Builds splined curve and defines key points
			//	@param distr: distribution to be splined
			//	@param base: moving window half-length
			//	@param summit[out]: returned X,Y coordinates of splined (smoothed) summit
			//	@returns key points: X-coord of highest point, X-coord of right middle hight point
			static fpair GetSplineKeyPoints(const dmap& distr, fraglen base, fpair& summit);

			// Calculates the best approximate distribution parameters for splined distribution
			//	@param distr: distribution to be splined
			//	@param dind: distribution index
			//	@param keypts[out]: returned key points
			void SetParamsForSpline(Distrib& distr, eDIndex dind, fpair& keypts);

			// Calculates the best approximate distribution parameters for interpolated distribution
			//	@param distr: distribution to be interpolated
			//	@param dind: distribution index
			//	@param keypts[out]: returned key points
			void SetParamsForInterpol(Distrib& distr, eDIndex dind, fpair& keypts);

			// Compares the real distribution with calculated one with given params, and sets PCC
			//	@param distr: compared real distribution
			//	@param dind: distribution index
			//	@param Mode: X-coordinate of summit
			//	@param full: if true then correlate from the beginning, otherwiase from summit
			//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
			void CalcPCC(const dmap& distr, eDIndex dind, int Mode, bool full = true);
		};

		// Returns combined distribution type by distribution index
		//static eDType GetCType(eDIndex dind) { return eDType(1 << dind); }

		// Returns distribution index by combined distribution type
		const static eDIndex GetDType(eDType dtype) { return eDIndex(RightOnePos(int(dtype))); }

		// Default constructor
		ADPs();

		float GetBestPCC() const { return _indADPs[0].PCC; }

		fpair GetBestParams() { Sort(); return _indADPs[0].Params; }

		// Calculates the best approximate distribution parameters for a specific type of distribution
		void CalcParams(eDType dtype, eDSmooth smode, Distrib& distr);

		// Sorts parameters and returns inner index of distribution with the highest PCC
		//	@returns inner index of distribution with the highest (best) PCC
		eDIndex GetBestIndex() { Sort(); return _indADPs[0].DIndex; }

		// Prints sorted distibutions params on a new line
		//	@param s: output stream
		void Print(dostream& s);

	private:
		// 'Indexed ADP': ADP supplied with inner index
		struct IndexedADP : public ADP
		{
			eDIndex	DIndex;
			bool	IsSet = false;	// true if this distribution approximation is required

			void Copy(const ADP& dp) { PCC = dp.PCC; Params = dp.Params; }

			// Prints AD parameters
			//	@param s: print stream
			//	@param maxPCC: maximum PCC to print relative PCC percentage
			void Print(dostream& s, float maxPCC) const;
		};

		std::array<IndexedADP, eDType::CNT>	_indADPs;
		bool _sorted = false;

		// Returns true if inner index is represented in combo cType
		static bool IsIndex(eDType cType, eDIndex dind) { return cType & (1 << dind); }

		// Returns true if AD parameters set in sorted instance
		bool IsSetInSorted(eDType dtype) const;

		// Returns number of ADP in sorted instance
		int SetSortedCount() const;

		// Returns ADP by combined distribution type
		//ADP& Params(eDType dtype) { return _indADPs[GetDType(dtype)]; }

		// Sorts in PCC descending order
		void Sort();

		// Sets distribution approximation computational requirement according to combo type
		//	@param dtype: required distribution types
		//	@returns: true if normal distribution is absent and added
		bool SetDTypes(eDType dtype);

		// Calculates and set approximate distribution parameters by index
		//	@param dind: distribution index
		//	@param adp: approximate distribution parameters
		void CopyParams(eDIndex dind, const ADP& adp) { _indADPs[dind].Copy(adp); }
	};

	// Returns specification string by specification type
	static const string Spec(eSpec s);

	eSpec _spec = eSpec::CLEAR;		// distribution specification
	ADPs	_indADPs;				// approximate distribution parameters for all type of distributions
	// these two fields are needed to print warnings after calculating the parameters
	fraglen	_base = FRAGLEN_MAX;	// moving window half-length of best spline
	//dpoint	_summit;				// X,Y coordinates of best splined (smoothed) summit
	fpair	_summit;				// X,Y coordinates of best splined (smoothed) summit
#ifdef MY_DEBUG
	mutable vector<dpoint> _spline;		// splining curve (container) to visualize splining
	mutable bool _fillSpline = true;	// true if fill splining curve (container)
	dostream* _s = NULL;				// print stream
#endif

	// Set moving window half-length of appropriate spline (estimated base)
	void EstimateSplineBase();

	// Prints warnings about poor quality distributions
	//	@param s: print stream
	void PrintWarning(dostream& s);

	// Prints original distribution as a set of <value>-<size> pairs
	//	@param s: print stream
	void PrintOriginal(dostream& s) const;
};
