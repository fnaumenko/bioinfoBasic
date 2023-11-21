/**********************************************************
Feature.h  2023 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 11/21/2023
***********************************************************/
#pragma once

#include "ChromData.h"
#include "DataInFile.h"

struct Featr : public Region
{
	float	Value;			// features's score

	Featr(const Region& rgn, float val = 0) : Value(val), Region(rgn) {}

	//Region& operator = (const Featr& f) { *this = f; }

#ifdef MY_DEBUG
	void Print() const { cout << Start << TAB << End << TAB << Value << LF; }
#endif	// MY_DEBUG
};

// 'Features' represents a collection of crhoms features
class Features : public Items<Featr>
{
	FBedInFile* _file = nullptr;		// valid only in constructor!
#ifdef _ISCHIP
	readlen	_minFtrLen;			// minimal length of feature
	float	_maxScore = 0;		// maximal feature score after reading
	bool	_uniScore = false;	// true if score is undefined in input data and set as 1
#elif defined _BIOCC
	// this is needed to get warning by calling Reads instead of Features (without -a option)
	bool	_narrowLenDistr = false;	// true if features length distribution is degenerate
#endif

	// Returns item's title
	//	@param pl: if true then print in plural form
	static const string& ItemTitle(bool pl = false) { return FT::ItemTitle(FT::eType::BED, pl); }

	// Returns a copy of Region
	//	@param it: item's iterator
	Region const Regn(cItemsIter it) const { return Region(*it); }

	// Adds chromosome to the instance
	//	@param cID: chromosome's ID
	//	@param cnt: count of chrom items
	void AddChrom(chrid cID, size_t cnt);

#ifdef _ISCHIP
	// Scales defined score through all features to the part of 1.
	void ScaleScores();
#endif

public:
#ifdef _ISCHIP
	// Creates new instance by bed-file name
	//	@param fName: file name
	//	@param cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@param joinOvrl: if true then join overlapping features, otherwise omit
	//	@param scoreInd: index of 'score' field
	//	@param bsLen: length of binding site: shorter features would be omitted
	//	@param prfName: true if file name should be printed unconditionally
	Features(const char* fName, ChromSizes& cSizes, bool joinOvrl,
		BYTE scoreInd, readlen bsLen, bool prfName)
		: _minFtrLen(bsLen), _uniScore(!scoreInd)
	{
		FBedInFile file(fName, &cSizes, scoreInd,
			joinOvrl ? UniBedInFile::eAction::JOIN : UniBedInFile::eAction::OMIT,
			eOInfo::LAC, prfName, true);
#else
	// Creates new instance by bed-file name
	//	@param fName: name of bed-file
	//	@param cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@param joinOvrl: if true then join overlapping features, otherwise omit
	//	@param prfName: true if file name should be printed unconditionally
	//	@param abortInvalid: true if invalid instance should abort excecution
	Features(const char* fName, ChromSizes & cSizes, bool joinOvrl,
		eOInfo oinfo, bool prfName, bool abortInvalid = true)
	{
		FBedInFile file(fName, &cSizes, 5,
			joinOvrl ? UniBedInFile::eAction::JOIN : UniBedInFile::eAction::OMIT,
			oinfo, prfName, abortInvalid);
#endif
		size_t estItemCnt = file.EstItemCount();
		if (estItemCnt) {
			ReserveItems(estItemCnt);
			_file = &file;
			file.Pass(*this);
			_file = nullptr;
		}
#ifdef _BIOCC
		_narrowLenDistr = file.NarrowLenDistr();
#endif
		//PrintEst(estItemCnt);
	}


	// Treats current item
	//	return: true if item is accepted
	bool operator()();

	// Closes current chrom, open next one
	//	@param cID: current chromosome's ID
	//	@param cLen: chromosome's  length
	//	@param cnt: current chromosome's items count
	//	@param nextcID: next chromosome's ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) { AddChrom(cID, cnt); }

	// Closes last chrom
	//	@param cID: last chromosome's ID
	//	@param cLen: chromosome's  length
	//	@param cnt: last chromosome's items count
	//	@param tCnt: total items count
	void operator()(chrid cID, chrlen cLen, size_t cnt, ULONG tCnt) { AddChrom(cID, cnt); }

	// Gets chromosome's feature by ID
	//	@param cID: chromosome's ID
	//	@param fInd: feature's index, or first feature by default
	//const Featr& Feature(chrid cID, chrlen fInd=0) const { return Item(cID, fInd); }

	// Gets chromosome's feature by iterator
	//	@param cit: chromosome's iterator
	//	@param fInd: feature's index, or first feature by default
	const Featr& Feature(cIter cit, chrlen fInd = 0) const { return Item(cit, fInd); }

	// Gets chromosome's feature region by iterator ???
	//	@param cit: chromosome's iterator
	//	@param fInd: feature's index, or first feature by default
	const Region& Regn(cIter cit, chrlen fInd = 0) const { return (const Region&)Item(cit, fInd); }

	// Gets the total length of all chromosome's features
	//	@param cit: chromosome's iterator
	chrlen FeaturesLength(cIter cit) const;

	// Gets chromosome's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@param cit: chromosome's iterator
	//	@param multiplier: 1 for numerics, 0 for letters
	//	@param fLen: average fragment length on which each feature will be expanded in puprose of calculation (float to minimize rounding error)
	chrlen EnrRegnLength(cIter cit, BYTE multiplier, float fLen) const {
		return (FeaturesLength(cit) + chrlen(2 * fLen * Data(cit).ItemsCount())) << multiplier;
	}

	// Gets chrom's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@param cID: chromosome's ID
	//	@param multiplier: 1 for numerics, 0 for nameds
	//	@param fLen: average fragment length on which each feature will be expanded in puprose of calculation (float to minimize rounding error)
	//	@return: chrom's total enriched regions length, or 0 if chrom is absent
	chrlen EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const;

	// Return min feature length
	chrlen GetMinFeatureLength() const;

	// Return min distance between features boundaries
	chrlen GetMinDistance() const;

	// Expands all features positions on the fixed length in both directions.
	// If extended feature starts from negative, or ends after chrom length, it is fitted.
	//	@param extLen: distance on which Start should be decreased, End should be increased, or inside out if it os negative
	//	@param cSizes: chrom sizes
	//	@param action: action for overlapping features
	//	@return: true if positions have been changed
	bool Extend(chrlen extLen, const ChromSizes & cSizes, UniBedInFile::eAction action);

	// Checks whether all features length exceed given length, throws exception otherwise.
	//	@param len: given control length
	//	@param lenDefinition: control length definition to print in exception message
	//	@param sender: exception sender to print in exception message
	void CheckFeaturesLength(chrlen len, const string & lenDefinition, const char* sender) const;

#ifdef _ISCHIP
	bool IsUniScore() const { return _uniScore; }
#else
	// Copies features coordinates to external DefRegions.
	void FillRegions(chrid cID, Regions & regn) const;
#endif	// _ISCHIP

#ifdef _BIOCC
	// Returns true if features length distribution is degenerate
	bool NarrowLenDistr() const { return _narrowLenDistr; }

	friend class JointedBeds;	// to access GetIter(chrid)
#endif

#ifdef MY_DEBUG
	void Print(size_t cnt = 0) const { Items::Print("features", cnt); }
#endif	// MY_DEBUG
};
