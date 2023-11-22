/**********************************************************
Distrib 2023 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 11/19/2023
-------------------------
Provides value (typically frequency) distribution functionality
***********************************************************/
#pragma once

#include "DataReader.h"
#include <array>

//#define MY_DEBUG

using dVal_t = size_t;
using fpair = pair<float, float>;

// 'Distrib' represents a value (typically fragment's/read's length frequency) distribution statistics
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

private:
	using dtype = int;	// consecutive distribution type: just to designate dist type, used as an index
	using spoint = pair<fraglen, dVal_t>;	// initial raw sequence point
	using dpoint = pair<fraglen, float>;	// distribution point 

	// Returns combined distribution type by consecutive distribution type
	static eCType GetCType(dtype type) { return eCType(1 << type); }

	// Returns consecutive distribution type by combined distribution type
	const static dtype GetDType(eCType ctype) { return RightOnePos(int(ctype)); }

	enum class eSpec {	// distribution specification
		CLEAR,		// normal quality;	exclusive
		SMOOTH,		// complementary
		MODUL,		// modulated;	complementary
		EVEN,		// exclusive
		CROP,		// cropped;	exclusive
		HCROP,		// heavily cropped;		exclusive
		SDEFECT,	// slightly defective;	exclusive
		DEFECT		// defective; exclusive
	};

	static const char* sTitle[];
	static const string sSpec[];
	static const string sParams;
	static const string sInaccurate;
	const fraglen smoothBase = 1;	// splining base for the smooth distribution
	static const int hRatio = 2;	// ratio of the summit height to height of the measuring point
	static const float lghRatio;	// log of ratio of the summit height to height of the measuring point

	// Keeps distribution params: PCC, mean(alpha), sigma(beta)
	struct DParams
	{
	private:
		static const float UndefPCC;
	public:
		float	PCC = 0;		// Pearson correlation coefficient
		fpair	Params{};		// mean(alpha), sigma(beta)

		bool operator >(const DParams& dp) const { return PCC > dp.PCC; }

		bool IsUndefPcc() const { return PCC == UndefPCC; };

		void SetUndefPcc() { PCC = UndefPCC; };
	};

	// 'AllDParams' represents a collection of restored distribution params for all type of distribution
	class AllDParams
	{
		// 'QualDParams' keeps restored distribution parameters
		struct QualDParams
		{
			eCType	Type;			//  combined distribution type
			const char* Title;		// distribution type title
			DParams dParams;		// PCC, mean(alpha), sigma(beta)

			// Sets combined type and title by consecutive type
			void SetTitle(dtype type) { Type = GetCType(type); Title = sTitle[type]; }

			// Returns true if distrib parameters set
			bool IsSet() const { return dParams.PCC != 0; }

			// Prints restored distr parameters
			//	@param s: print stream
			//	@param maxPCC: masimum PCC to print relative PCC percentage
			void Print(dostream& s, float maxPCC) const;
		};

		array<QualDParams, eCType::CNT>	_allParams;
		bool _sorted = false;

		// Returns true if distribution parameters set in sorted instance
		bool IsSetInSorted(eCType ctype) const;

		// Returns number of distribution parameters set in sorted instance
		int SetCntInSorted() const;

		// Returns DParams by combined distribution type
		DParams& Params(eCType ctype) { return _allParams[GetDType(ctype)].dParams; }

		// Sorts in PCC descending order
		void Sort();

	public:
		// Default constructor
		AllDParams();

		// Set distribution parameters by type
		//	@param type: consecutive distribution type
		//	@param dp: PCC, mean(alpha) & sigma(beta)
		void SetParams(dtype type, const DParams& dp) { _allParams[type].dParams = dp; }

		// Clear normal distribution if its PCC is less then lognorm PCC by the threshold
		void ClearNormDistBelowThreshold(float thresh) {
			if (Params(eCType::LNORM).PCC / Params(eCType::NORM).PCC > thresh)
				Params(eCType::NORM).PCC = 0;
		}

		// Returns parameters of distribution with the highest (best) PCC
		//	@param dParams: returned PCC, mean(alpha) & sigma(beta)
		//	@return: consecutive distribution type with the highest (best) PCC
		dtype GetBestParams(DParams& dParams);

		// Prints sorted distibutions params
		//	@param s: print stream
		void Print(dostream& s);
	};

	// Returns specification string by specification type
	static const string Spec(eSpec s) { return "Distribution " + sSpec[int(s)]; }

	// Returns true if consecutive type is represented in combo cType
	static bool IsType(eCType cType, dtype type) { return cType & (1 << type); }

	// Returns true if combo type is represented in combo cType
	static bool IsType(eCType cType, eCType type) { return cType & type; }

	// Calls distribution parameters by consecutive distribution type
	//	@param keypts: key points: X-coord of highest point, X-coord of right middle hight point
	//	@param params: returned mean(alpha) & sigma(beta)
	//	Defined as a member of the class only to use the private short name lghRatio
	static void (*SetParams[])(const fpair& keypts, fpair& params);

	eSpec _spec = eSpec::CLEAR;		// distribution specification
	AllDParams	_allParams;			// distributions parameters
#ifdef MY_DEBUG
	mutable vector<dpoint> _spline;		// splining curve (container) to visualize splining
	mutable bool _fillSpline = true;	// true if fill splining curve (container)
	dostream* _s = NULL;				// print stream
#endif

	// Returns estimated moving window half-length ("base")
	//	@return: estimated base, or 0 in case of degenerate distribution
	fraglen GetBase();

	// Defines key points
	//	@param base: moving window half-length
	//	@param summit: returned X,Y coordinates of spliced (smoothed) summit
	//	@return: key points: X-coord of highest point, X-coord of right middle hight point
	fpair GetKeyPoints(fraglen base, dpoint& summit) const;

	// Compares this sequence with calculated one by given mean & sigma, and returns PCC
	//	@param type: consecutive distribution type
	//	@param dParams: returned PCC, input mean(alpha) & sigma(beta)
	//	@param Mode: X-coordinate of summit
	//	@param full: if true then correlate from the beginning, otherwiase from summit
	//	calculated on the basis of the "start of the sequence" – "the first value less than 0.1% of the maximum".
	void CalcPCC(dtype type, DParams& dParams, fraglen Mode, bool full = true) const;

	// Calculates distribution parameters
	//	@param type: consecutive distribution type
	//	@param keyPts: key points: X-coord of highest point, X-coord of right middle hight point
	//	@param dParams: returned PCC, mean(alpha) & sigma(beta)
	//	@param Mode: X-coordinate of summit
	void SetPCC(dtype type, const fpair& keypts, DParams& dParams, fraglen Mode) const;

	// Calculates and print called distribution parameters
	//	@param type: consecutive distribution type
	//	@param base: moving window half-length
	//	@param summit: returned X,Y coordinates of best spliced (smoothed) summit
	void CallParams(dtype type, fraglen base, dpoint& summit);

	// Prints original distribution features
	//	@param s: print stream
	//	@param base: moving window half-length
	//	@param summit: X,Y coordinates of spliced (smoothed) summit
	void PrintTraits(dostream& s, fraglen base, const dpoint& summit);

	// Prints original distribution as a set of <value>-<size> pairs
	//	@param s: print stream
	void PrintSeq(dostream& s) const;

public:
	// Default constructor
	Distrib() {}

	// Constructor by ready distribution file
	//	@param fname: name of ready distribution file
	Distrib(const char* fname);

	// Returns true if distribution has not enough size
	//bool IsDegenerate() const { return size() < 5; }

	// Adds value to the instance
	void AddVal(fraglen val) { (*this)[val]++; }

	// Calculate and print distribution
	//	@param s: print stream
	//	@param type: combined type of distribution
	//	@param prDistr: if true then print distribution additionally
	void Print(dostream& s, eCType type, bool prDistr = true);
};
