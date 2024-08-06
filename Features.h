/**********************************************************
Feature.h
BED feature and features collection
Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 08/06/2024
***********************************************************/
#pragma once

#include "ChromData.h"
#include "DataReader.h"

// 'ItemIndices' representes a range of chromosome's indices
struct ItemIndices
{
	size_t	FirstInd;	// first index in items container
	size_t	LastInd;	// last index in items container

	ItemIndices(size_t first = 0, size_t last = 1) : FirstInd(first), LastInd(last - 1) {}

	// Returns count of items
	size_t ItemsCount() const { return LastInd - FirstInd + 1; }
};

template <typename I>
class Items : public Chroms<ItemIndices>
{
public:
	typedef typename vector<I>::const_iterator cItemsIter;

protected:
	typedef typename vector<I>::iterator ItemsIter;

	vector<I> _items;

	// Applies a function to each item of the specified chromosome
	//	@param cit: chromosome's iterator
	//	@param fn: function accepted item iterator
	void DoForChrItems(cIter cit, function<void(cItemsIter)> fn) const {
		const auto itEnd = ItemsEnd(Data(cit));
		for (auto it = ItemsBegin(Data(cit)); it != itEnd; it++)
			fn(it);
	}

	// Applies a function to each item of the whole collection
	//	@param fn: function accepted item iterator
	void DoForItems(function<void(cItemsIter)> fn) const {
		for (const auto& c : Container()) {
			const auto itEnd = ItemsEnd(c.second.Data);
			for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
				fn(it);
		}
	}

	// Initializes size of positions container.
	void ReserveItems(size_t size) { _items.reserve(size); }

	// Returns item by chrom's iterator
	//	@param cit: chromosome's iterator
	//	@param iInd: index of item
	const I& Item(cIter cit, chrlen iInd) const { return _items[Data(cit).FirstInd + iInd]; }

	// Returns item by chrom's ID
	//	@param cID: chromosome's ID
	//	@param iInd: index of item
	//const I& Item(chrid cID, chrlen iInd) const { return Item(GetIter(cID), iInd); }

#ifdef MY_DEBUG
	void PrintItem(chrlen itemInd) const { _items[itemInd].Print(); }
#endif

public:
	// === items count

	// Gets total count of items
	size_t ItemsCount() const { return _items.size(); }

	// Gets count of items for chrom by chrom's data
	//	@param data: chromosome's data
	size_t ItemsCount(const ItemIndices& data) const { return data.ItemsCount(); }

	// Gets count of items for chrom by chrom's iterator
	//	@param cit: chromosome's iterator
	size_t ItemsCount(cIter cit) const { return ItemsCount(Data(cit)); }

	// Gets count of items for chrom by chrom's ID
	//	@param cID: chromosome's ID
	size_t ItemsCount(chrid cID) const { return ItemsCount(GetIter(cID)); }

	// === iterator

	// Returns a constant iterator referring to the first item of specified item indices
	//	@param data: item indices
	cItemsIter ItemsBegin(const ItemIndices& data) const { return _items.begin() + data.FirstInd; }

	// Returns a constant iterator referring to the past-the-end item of specified item indices
	//	@param data: item indices
	cItemsIter ItemsEnd(const ItemIndices& data) const { return _items.begin() + data.LastInd + 1; }

	// Returns a constant iterator referring to the first item of specified chromosome's iterator
	//	@param cit: chromosome's iterator
	cItemsIter ItemsBegin(cIter cit) const { return ItemsBegin(Data(cit)); }

	// Returns a constant iterator referring to the past-the-end item of specified chromosome's iterator
	//	@param cit: chromosome's iterator
	cItemsIter ItemsEnd(cIter cit) const { return ItemsEnd(Data(cit)); }

	// Returns a constant iterator referring to the first item of specified chromosome's ID
	//	@param cID: chromosome's ID
	//cItemsIter ItemsBegin(chrid cID) const { return ItemsBegin(GetIter(cID)); }	// ineffective

	// Returns a constant iterator referring to the past-the-end item of specified chromosome's ID
	//	@param cID: chromosome's ID
	//cItemsIter ItemsEnd(chrid cID) const { return ItemsEnd(GetIter(cID)); }		// ineffective


	// Returns an iterator referring to the first item of specified item indices
	//	@data: item indices
	ItemsIter ItemsBegin(ItemIndices& data) { return _items.begin() + data.FirstInd; }

	// Returns an iterator referring to the past-the-end item of specified item indices
	//	@param data: item indices
	ItemsIter ItemsEnd(ItemIndices& data) { return _items.begin() + data.LastInd + 1; }

	// Returns an iterator referring to the first item of specified chromosome's iterator
	//	@param cit: chromosome's iterator
	//ItemsIter ItemsBegin(Iter cit) { return ItemsBegin(Data(cit)); }

	// Returns an iterator referring to the past-the-end item of specified chromosome's iterator
	//	@param cit: chromosome's iterator
	//ItemsIter ItemsEnd(Iter cit) { return ItemsEnd(Data(cit)); }

	// === print

	void PrintEst(ULONG estCnt) const { cout << " est/fact: " << double(estCnt) / ItemsCount() << LF; }

#ifdef _ISCHIP
	// Prints items name and count, adding chrom name if the instance holds only one chrom
	//	@param ftype: file type
	//	@param prLF: if true then print line feed at the end
	void PrintItemCount(FT::eType ftype, bool prLF = true) const
	{
		size_t iCnt = ItemsCount();
		dout << iCnt << SPACE << FT::ItemTitle(ftype, iCnt > 1);
		if (ChromCount() == 1)		dout << " per " << Chrom::TitleName(CID(cBegin()));
		if (prLF)	dout << LF;
	}
#endif

#ifdef MY_DEBUG
	// Prints collection
	//	@param title: item title
	//	@param prICnt: number of printed items per each chrom, or 0 if all
	void Print(const char* title, size_t prICnt = 0) const
	{
		cout << LF << title << SepCl;
		if (prICnt)	cout << "first " << prICnt << " per each chrom\n";
		else		cout << ItemsCount() << LF;
		for (const auto& inds : Container()) {
			const string& chr = Chrom::AbbrName(inds.first);
			const auto& data = inds.second.Data;
			const size_t lim = prICnt ? prICnt + data.FirstInd : UINT_MAX;
			for (chrlen i = chrlen(data.FirstInd); i <= chrlen(data.LastInd); i++) {
				if (i >= lim)	break;
				cout << chr << TAB;
				PrintItem(i);
			}
		}
	}
#endif	// MY_DEBUG
};


struct Featr : public Region
{
	float	Value;			// features's score

	Featr(const Region& rgn, float val = 0) : Value(val), Region(rgn) {}

	//Region& operator = (const Featr& f) { *this = f; }

#ifdef MY_DEBUG
	void Print() const { cout << Start << TAB << End << TAB << Value << LF; }
#endif
};

// 'Features' represents a collection of crhoms features
class Features : public Items<Featr>
{
	FBedReader* _file = nullptr;	// valid only in constructor
#ifdef _FEATR_SCORE
	float	_maxScore = 0;			// maximal feature score after reading
	bool	_uniScore = false;		// true if score is undefined in input data and set as 1
#endif
#ifdef _ISCHIP
	readlen	_minFtrLen;				// minimal length of feature
#elif defined _BIOCC
	// this is needed to get warning by calling Reads instead of Features (without -a option)
	bool	_narrowLenDistr = false;	// true if features length distribution is degenerate
#endif

	// Returns item's title
	//	@param pl: if true then print in plural form
	static const string& ItemTitle(bool pl = false) { return FT::ItemTitle(FT::BED, pl); }

	// Returns a copy of Region
	//	@param it: item's iterator
	Region const Regn(cItemsIter it) const { return Region(*it); }

	// Initializes new instance by bed-file name
	//	@param fName: name of bed-file
	//	@param cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@param scoreInd: index of 'score' field
	//	@param joinOvrl: if true then join overlapping features, otherwise omit
	//	@param oinfo: outputted info
	//	@param abortInvalid: true if invalid instance should abort excecution
	void Init(
		const char* fName,
		ChromSizes* cSizes,
		BYTE scoreInd,
		bool joinOvrl,
		eOInfo oinfo,
		bool abortInvalid
	);

	// Adds chromosome to the instance
	//	@param cID: chromosome's ID
	//	@param cnt: count of chrom items
	void AddChrom(chrid cID, size_t cnt);

	// Scales defined score through all features to the part of 1.
	//void ScaleScores();

public:
#ifdef _ISCHIP
	// Creates new instance by bed-file name
	//	@param fName: file name
	//	@param cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@param scoreInd: index of 'score' field
	//	@param joinOvrl: if true then join overlapping features, otherwise omit
	//	@param bsLen: length of binding site: shorter features would be omitted
	Features(const char* fName, ChromSizes& cSizes, BYTE scoreInd, bool joinOvrl, readlen bsLen)
		: _minFtrLen(bsLen), _uniScore(!scoreInd)
	{
		Init(fName, &cSizes, scoreInd, joinOvrl, eOInfo::LAC, true);
	}
#else
	// Creates new instance by bed-file name
	//	@param fName: name of bed-file
	//	@param cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@param joinOvrl: if true then join overlapping features, otherwise omit
	//	@param oinfo: outputted info
	//	@param abortInvalid: true if invalid instance should abort excecution
	Features(const char* fName, ChromSizes* cSizes, bool joinOvrl,
		eOInfo oinfo, bool abortInvalid = true)
	{
		Init(fName, cSizes, 5, joinOvrl, oinfo, abortInvalid);
	}
#endif

	// Treats current item
	//	@returns: true if item is accepted
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
	void operator()(chrid cID, chrlen cLen, size_t cnt, size_t tCnt) { AddChrom(cID, cnt); }

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
	//	@returns: chrom's total enriched regions length, or 0 if chrom is absent
	chrlen EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const;

#ifdef _ISCHIP
	bool IsUniScore() const { return _uniScore; }

	// Return minimal feature length
	chrlen GetMinFeatureLength() const;
#endif	// _ISCHIP

#ifdef _BIOCC
	// Returns true if features length distribution is degenerate
	bool NarrowLenDistr() const { return _narrowLenDistr; }

	friend class JointedBeds;	// to access GetIter(chrid)
#endif

	// Return min distance between features boundaries
	chrlen GetMinDistance() const;

	// Increases the size of each feature in both directions.
	// If expanded feature starts from negative, or ends after chrom length, it is fitted.
	//	@param expLen: value on which Start should be decreased, End should be increased
	//	@param cSizes: chromosome sizes for chromosome's length control or NULL if no control
	//	@param action: action for overlapping features
	//	@returns: true if the expansion was completed successfully
	bool Expand(chrlen expLen, const ChromSizes* cSizes, UniBedReader::eAction action);

	// Checks whether all features length exceed given length, throws exception otherwise.
	//	@param len: given control length
	//	@param lenDefinition: control length definition to print in exception message
	//	@param sender: exception sender to print in exception message
	void CheckFeaturesLength(chrlen len, const string & lenDefinition, const char* sender) const;

	// Copies features coordinates to external DefRegions.
	//void FillRegions(chrid cID, Regions& regn) const;

#ifdef MY_DEBUG
	void Print(size_t cnt = 0) const { Items::Print("features", cnt); }
#endif	// MY_DEBUG
};
