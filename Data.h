/**********************************************************
Data.h (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 11/18/2023
-------------------------
Provides common data functionality
***********************************************************/
#pragma once

#include "DataInFile.h"
#include <array>
#include <algorithm>    // std::sort

#define	CID(it)	(it)->first

//static const string range_out_msg = "ChromMap[]: invalid key ";

// 'ChromMap'
template <typename T> class ChromMap
{
public:
	typedef map<chrid, T> chrMap;
	typedef typename chrMap::iterator Iter;			// iterator
	typedef typename chrMap::const_iterator cIter;	// constant iterator

	// Returns a random-access constant iterator to the first element in the container
	cIter cBegin() const { return _cMap.begin(); }
	// Returns the past-the-end constant iterator.
	cIter cEnd()	  const { return _cMap.end(); }
	// Returns a random-access iterator to the first element in the container
	Iter Begin() { return _cMap.begin(); }
	// Returns the past-the-end iterator.
	Iter End() { return _cMap.end(); }

	// Copies entry
	void Assign(const ChromMap& map) { _cMap = map._cMap; }

	// Returns constant reference to the item at its cID
	const T& At(chrid cID) const { return _cMap.at(cID); }

	// Returns reference to the item at its cID
	T& At(chrid cID) { return _cMap.at(cID); }

	const T& operator[] (chrid cID) const { return _cMap(cID); }

	T& operator[] (chrid cID) { return _cMap[cID]; }

	// Searches the container for a cID and returns an iterator to the element if found,
	// otherwise it returns an iterator to end (the element past the end of the container)
	Iter GetIter(chrid cID) { return _cMap.find(cID); }

	// Searches the container for a cID and returns a constant iterator to the element if found,
	// otherwise it returns an iterator to cEnd (the element past the end of the container)
	const cIter GetIter(chrid cID) const { return _cMap.find(cID); }

	// Adds value type to the collection without checking cID
	void AddVal(chrid cID, const T& val) { _cMap[cID] = val; }
	void AddVal(chrid cID, const T&& val) { _cMap[cID] = val; }

	//// Adds value type to the collection without checking cID
	//void EmplaceVal(chrid cID, BYTE val) { _cMap.emplace(cID, T(val)); }
	////void EmplaceVal(T&&... args) { _cMap.emplace(cID, args); }


	// Removes from the container an element cID
	void Erase(chrid cID) { _cMap.erase(cID); }

	// Clear content
	void Clear() { _cMap.clear(); }

	// Returns true if element with cID exists in the container, and false otherwise.
	bool FindItem(chrid cID) const { return _cMap.count(cID) > 0; }

protected:
	// Returns count of elements.
	size_t Count() const { return _cMap.size(); }

	// Adds class type to the collection without checking cID.
	// Avoids unnecessery copy constructor call
	//	return: class type collection reference
	T& AddElem(chrid cID, const T& val) { return _cMap[cID] = val; }

	// Adds empty class type to the collection without checking cID
	//	return: class type collection reference
	T& AddEmptyElem(chrid cID) { return AddElem(cID, T()); }

	const chrMap& Container() const { return _cMap; }

	chrMap& Container() { return _cMap; }

private:
	chrMap _cMap;
};

// 'ChromData' implements sign whether chrom is involved in processing, and chrom's data itself
template <typename T> struct ChromData
{
	bool Treated = true;	// true if chrom is involved in processing
	T	 Data;				// chrom's data

	ChromData() : Data(T()) {}
	ChromData(const T& data) : Data(data) {}
};

// Basic class for chromosomes collection; keyword 'abstract' doesn't compiled in gcc
template <typename T>
class Chroms : public ChromMap<ChromData<T> >
{
public:
	// Returns true if chromosome by iterator should be treated
	//bool IsTreated(typename ChromMap<ChromData<T> >& data) const { return data.second.Treated; }

	// Returns true if chromosome by iterator should be treated
	bool IsTreated(typename ChromMap<ChromData<T> >::cIter it) const { return it->second.Treated; }

	// Returns true if chromosome by ID should be treated
	//bool IsTreated(chrid cID) const { return IsTreated(this->GetIter(cID));	}

	const T& Data(typename ChromMap<ChromData<T> >::cIter it) const { return it->second.Data; }

	T& Data(typename ChromMap<ChromData<T> >::Iter it) { return it->second.Data; }

	T& Data(chrid cID) { return this->At(cID).Data; }

	// Returns count of chromosomes.
	chrid ChromCount() const { return chrid(this->Count()); }	// 'this' required in gcc

	 // Gets count of treated chroms
	chrid TreatedCount() const
	{
		chrid res = 0;
		for (auto it = this->cBegin(); it != this->cEnd(); res += IsTreated(it++));
		return res;
	}

	// Returns true if chromosome cID exists in the container, and false otherwise.
	bool FindChrom(chrid cID) const { return this->FindItem(cID); }

	// Adds value type to the collection without checking cID
	void AddValue(chrid cid, const T& val) { this->AddVal(cid, ChromData<T>(val)); }

#ifdef	_READDENS
	// Sets common chromosomes as 'Treated' in both of this instance and given Chroms.
	// Objects in both collections have to have second field as Treated.
	//	@obj: compared Chroms object
	//	@printWarn: if true print warning - uncommon chromosomes
	//	@throwExcept: if true then throw exception if no common chroms finded, otherwise print warning
	//	return: count of common chromosomes
	chrid	SetCommonChroms(Chroms<T>& obj, bool printWarn, bool throwExcept)
	{
		typename ChromMap<ChromData<T> >::Iter it;
		chrid commCnt = 0;

		// set treated chroms in this instance
		for (it = this->Begin(); it != this->End(); it++)	// 'this' required in gcc
			if (it->second.Treated = obj.FindChrom(CID(it)))
				commCnt++;
			else if (printWarn)
				Err(Chrom::Absent(CID(it), "second file")).Warning();
		// set false treated chroms in a compared object
		for (it = obj.Begin(); it != obj.End(); it++)
			if (!FindChrom(CID(it))) {
				it->second.Treated = false;
				if (printWarn)
					Err(Chrom::Absent(CID(it), "first file")).Warning();
			}
		if (!commCnt)	Err("no common " + Chrom::Title(true)).Throw(throwExcept);
		return commCnt;
	}
#endif	// _READDENS
};

// 'ItemIndices' representes a range of chromosome's features/reads indexes,
struct ItemIndices
{
	chrlen	FirstInd;	// first index in items container
	chrlen	LastInd;	// last index in items container

	ItemIndices(chrlen first = 0, chrlen last = 1) : FirstInd(first), LastInd(last - 1) {}

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

	vector<I> _items;	// vector of bed items

	// Applies function fn to each of the item for chrom defined by cit
	void ForChrItems(cIter cit, function<void(cItemsIter)> fn) const {
		const auto itEnd = ItemsEnd(Data(cit));
		for (auto it = ItemsBegin(Data(cit)); it != itEnd; it++)
			fn(it);
	}

	// Applies function fn to each of the item for chrom defined by cit
	//void ForChrItems(const value_type& c, function<void(cItemsIter)> fn) const {
	//	//const ItemIndices& data = 
	//	const auto itEnd = ItemsEnd(c.Data.
	//		//Data(cit));
	//	for (auto it = ItemsBegin(Data(cit)); it != itEnd; it++)
	//		fn(it);
	//}

	// Applies function fn to each of the item for all chroms
	void ForAllItems(function<void(cItemsIter)> fn) const {
		for (const auto& c : Container()) {
			const auto itEnd = ItemsEnd(c.second.Data);
			for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
				fn(it);
		}
	}

	// Initializes size of positions container.
	void ReserveItems(size_t size) { _items.reserve(size); }

	// Gets item.
	//	@it: chromosome's iterator
	//	@iInd: index of item
	const I& Item(cIter it, chrlen iInd) const { return _items[Data(it).FirstInd + iInd]; }

	// Returns item.
	//	@cInd: index of chromosome
	//	@iInd: index of item
	//const I& Item(chrid cID, chrlen iInd) const { return Item(GetIter(cID), iInd); }

#ifdef _DEBUG
	void PrintItem(chrlen itemInd) const { _items[itemInd].Print(); }
#endif

public:
	// Gets total count of items
	size_t ItemsCount() const { return _items.size(); }

	// Gets count of items for chrom
	//	@data: item indexes
	size_t ItemsCount(const ItemIndices& data) const { return data.ItemsCount(); }

	// Gets count of items for chrom
	//	@cit: chrom's iterator
	size_t ItemsCount(cIter cit) const { return ItemsCount(Data(cit)); }

	// Gets count of items for chrom; used in isChIP
	//	@cit: chrom's ID
	size_t ItemsCount(chrid cID) const { return ItemsCount(GetIter(cID)); }

	void PrintEst(ULONG estCnt) const { cout << " est/fact: " << double(estCnt) / ItemsCount() << LF; }

#ifdef _ISCHIP
	// Prints items name and count, adding chrom name if the instance holds only one chrom
	//	@ftype: file type
	//	@prLF: if true then print line feed
	void PrintItemCount(FT::eType ftype, bool prLF = true) const
	{
		size_t iCnt = ItemsCount();
		dout << iCnt << SPACE << FT::ItemTitle(ftype, iCnt > 1);
		if (ChromCount() == 1)		dout << " per " << Chrom::TitleName(CID(cBegin()));
		if (prLF)	dout << LF;
	}
#endif
	// Returns a constan iterator referring to the first item of specified chrom
	//	@data: item indexes
	cItemsIter ItemsBegin(const ItemIndices& data) const { return _items.begin() + data.FirstInd; }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@data: item indexes
	cItemsIter ItemsEnd(const ItemIndices& data) const { return _items.begin() + data.LastInd + 1; }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@cit: chromosome's constant iterator
	cItemsIter ItemsBegin(cIter cit) const { return ItemsBegin(Data(cit)); }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@cID: chromosome's ID
	cItemsIter ItemsBegin(chrid cID) const { return ItemsBegin(GetIter(cID)); }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome'sconstant  iterator
	cItemsIter ItemsEnd(cIter cit) const { return ItemsEnd(Data(cit)); }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@cID: chromosome's ID
	cItemsIter ItemsEnd(chrid cID) const { return ItemsEnd(GetIter(cID)); }

	// Returns a constan iterator referring to the first item of specified chrom
	//	@data: item indexes
	ItemsIter ItemsBegin(ItemIndices& data) { return _items.begin() + data.FirstInd; }

	// Returns a constant iterator referring to the past-the-end item of specified chrom
	//	@data: item indexes
	ItemsIter ItemsEnd(ItemIndices& data) { return _items.begin() + data.LastInd + 1; }

	// Returns an iterator referring to the first item of specified chrom
	//	@cit: chromosome's iterator
	ItemsIter ItemsBegin(Iter cit) { return ItemsBegin(Data(cit)); }

	// Returns an iterator referring to the past-the-end item of specified chrom
	//	@cit: chromosome's iterator
	ItemsIter ItemsEnd(Iter cit) { return ItemsEnd(Data(cit)); }

#ifdef _DEBUG
	// Prints collection
	//	@title: item title
	//	@prICnt: number of printed items per each chrom, or 0 if all
	void Print(const char* title, size_t prICnt = 0) const
	{
		cout << LF << title << ": ";
		if (prICnt)	cout << "first " << prICnt << " per each chrom\n";
		else		cout << ItemsCount() << LF;
		for (const auto& inds : Container()) {
			const string& chr = Chrom::AbbrName(inds.first);
			const auto& data = inds.second.Data;
			const size_t lim = prICnt ? prICnt + data.FirstInd : UINT_MAX;
			for (chrlen i = data.FirstInd; i <= data.LastInd; i++) {
				if (i >= lim)	break;
				cout << chr << TAB;
				PrintItem(i);
			}
		}
	}
#endif	// _DEBUG
};

#ifdef _FEATURES

struct Featr : public Region
{
	float	Value;			// features's score

	Featr(const Region& rgn, float val = 0) : Value(val), Region(rgn) {}

	//Region& operator = (const Featr& f) { *this = f; }
#ifdef _DEBUG
	void Print() const { cout << Start << TAB << End << TAB << Value << LF; }
#endif	// _DEBUG
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

	// Gets item's title
	// Obj abstract method implementation.
	//	@pl: true if plural form
	const string& ItemTitle(bool pl = false) const { return FT::ItemTitle(FT::eType::BED, pl); }

	// Gets a copy of Region by item's iterator.
	// Abstract BaseItems<> method implementation.
	Region const Regn(cItemsIter it) const { return Region(*it); }

	// Adds chrom to the instance
	//	@cID: chrom
	//	@cnt: count of chrom items
	void AddChrom(chrid cID, size_t cnt);

#ifdef _ISCHIP
	// Scales defined score through all features to the part of 1.
	void ScaleScores();
#endif

public:
#ifdef _ISCHIP
	// Creates new instance by bed-file name
	//	@fName: file name
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@joinOvrl: if true then join overlapping features, otherwise omit
	//	@scoreInd: index of 'score' field
	//	@bsLen: length of binding site: shorter features would be omitted
	//	@prfName: true if file name should be printed unconditionally
	Features(const char* fName, ChromSizes& cSizes, bool joinOvrl,
		BYTE scoreInd, readlen bsLen, bool prfName)
		: _minFtrLen(bsLen), _uniScore(!scoreInd)
	{
		FBedInFile file(fName, &cSizes, scoreInd,
			joinOvrl ? UniBedInFile::eAction::JOIN : UniBedInFile::eAction::OMIT,
			eOInfo::LAC, prfName, true);
#else
	// Creates new instance by bed-file name
	//	@fName: name of bed-file
	//	@cSizes: chrom sizes to control the chrom length exceedeng, or NULL if no control
	//	@joinOvrl: if true then join overlapping features, otherwise omit
	//	@prfName: true if file name should be printed unconditionally
	//	@abortInvalid: true if invalid instance should abort excecution
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


	// treats current item
	//	return: true if item is accepted
	bool operator()();

	// Closes current chrom, open next one
	//	@cID: current chrom ID
	//	@cLen: chrom length
	//	@cnt: current chrom items count
	//	@nextcID: next chrom ID
	void operator()(chrid cID, chrlen cLen, size_t cnt, chrid nextcID) { AddChrom(cID, cnt); }

	// Closes last chrom
	//	@cID: last chrom ID
	//	@cLen: chrom length
	//	@cnt: last chrom items count
	//	@tCnt: total items count
	void operator()(chrid cID, chrlen cLen, size_t cnt, ULONG tCnt) { AddChrom(cID, cnt); }

	// Gets chromosome's feature by ID
	//	@cID: chromosome's ID
	//	@fInd: feature's index, or first feature by default
	//const Featr& Feature(chrid cID, chrlen fInd=0) const { return Item(cID, fInd); }

	// Gets chromosome's feature by iterator
	//	@it: chromosome's iterator
	//	@fInd: feature's index, or first feature by default
	const Featr& Feature(cIter it, chrlen fInd = 0) const { return Item(it, fInd); }

	// Gets chromosome's feature by iterator
	//	@it: chromosome's iterator
	//	@fInd: feature's index, or first feature by default
	const Region& Regn(cIter it, chrlen fInd = 0) const { return (const Region&)Item(it, fInd); }

	// Gets the sum length of all chromosome's features
	//	@it: chromosome's iterator
	chrlen FeaturesLength(cIter it) const;

	// Gets chromosome's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@it: chromosome's iterator
	//	@multiplier: 1 for numerics, 0 for letters
	//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
	//	(float to minimize rounding error)
	chrlen EnrRegnLength(cIter it, BYTE multiplier, float fLen) const {
		return (FeaturesLength(it) + chrlen(2 * fLen * Data(it).ItemsCount())) << multiplier;
	}

	// Gets chrom's total enriched regions length:
	// a double length for numeric chromosomes or a single for named.
	//	@cID: chromosome's ID
	//	@multiplier: 1 for numerics, 0 for nameds
	//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
	//	(float to minimize rounding error)
	//	return: chrom's total enriched regions length, or 0 if chrom is absent
	chrlen EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const;

	// Return min feature length
	chrlen GetMinFeatureLength() const;

	// Return min distance between features boundaries
	chrlen GetMinDistance() const;

	// Expands all features positions on the fixed length in both directions.
	// If extended feature starts from negative, or ends after chrom length, it is fitted.
	//	@extLen: distance on which Start should be decreased, End should be increased,
	//	or inside out if it os negative
	//	@cSizes: chrom sizes
	//	@action: action for overlapping features
	//	return: true if positions have been changed
	bool Extend(chrlen extLen, const ChromSizes & cSizes, UniBedInFile::eAction action);

	// Checks whether all features length exceed given length, throws exception otherwise.
	//	@len: given control length
	//	@lenDefinition: control length definition to print in exception message
	//	@sender: exception sender to print in exception message
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
#ifdef _DEBUG
	void Print(size_t cnt = 0) const { Items::Print("features", cnt); }
#endif	// _DEBUG
	};

#endif	// _FEATURES

// 'ChromSize' represents real and defined effective chrom lengths
struct ChromSize
{
	chrlen Real;				// real (actual) chrom length
#ifdef _ISCHIP
	mutable chrlen Defined = 0;	// defined effective chrom length;
								// 'effective' means double length for autosomes, single one for somatic

	// Sets chrom's effective (treated) real length as defined
	chrlen SetEffDefined(bool autosome) const { return bool(Defined = (Real << int(autosome))); }
#endif

	ChromSize(chrlen size = 0) : Real(size) {}
};

// 'ChromSizes' represented chrom sizes with additional file system binding attributes
// Holds path to reference genome and to service files
class ChromSizes : public Chroms<ChromSize>
{
	string	_ext;			// FA files real extention; if empty then instance is initialized by service dir
	string	_gPath;			// ref genome path
	string	_sPath;			// service path
	mutable genlen _gsize;	// size of whole genome


	// Returns length of common prefix before abbr chrom name of all file names
	//	@fName: full file name
	//	@extLen: length of file name's extention
	//	return: length of common prefix or -1 if there is no abbreviation chrom name in fName
	static int	CommonPrefixLength(const string& fName, BYTE extLen);

	// Initializes chrom sizes from file
	void Read(const string& fName);

	// Saves chrom sizes to file
	//	@fName: full file name
	void Write(const string& fName) const;

	// Fills external vector by chrom IDs relevant to file's names found in given directory.
	//	@cIDs: filling vector of chrom's IDs
	//	@gName: path to reference genome
	//	return: count of filled chrom's IDs
	chrid GetChromIDs(vector<chrid>& cIDs, const string& gName);

	// Initializes the paths
	//	@gPath: reference genome directory
	//	@sPath: service directory
	//	@prMsg: true if print message about service fodler and chrom.sizes generation
	void SetPath(const string& gPath, const char* sPath, bool prMsg);

	void SetTreatedChrom(chrid cID);

	// returns true if service path is defined
	bool IsServAvail() const { return _sPath.size(); }

protected:
	chrlen Length(cIter it) const { return Data(it).Real; }

public:
	const string& RefExt() const { return _ext; }

	// Creates and initializes an instance.
	//	@param gName: reference genome directory or chrom.sizes file
	//	@param customChrOpt: id of 'custom chrom' option
	//	@param prMsg: if true then print message about service fodler and chrom.sizes generation
	//	@param sPath: service directory
	//	@param checkGRef: if true then check if @gName is a ref genome dir; used in isChIP
	ChromSizes(const char* gName, BYTE customChrOpt, bool prMsg, const char* sPath = NULL, bool checkGRef = false);

	ChromSizes() { _ext = _gPath = _sPath = strEmpty; }

	// Initializes ChromSizes by SAM header
	void Init(const string& headerSAM);

	bool IsFilled() const { return Count(); }

	// Return true if chrom.sizes are defined explicitly, by user
	//bool IsExplicit() const { return !_gPath.length(); }		// false if path == strEmpty

	// Returns reference directory
	const string& RefPath() const { return _gPath; }

	// Returns service directory
	const string& ServPath() const { return _sPath; }

	// Returns service directory
	const bool IsServAsRef() const { return _gPath == _sPath; }

	// Returns full ref chrom name by chrom ID 
	const string RefName(chrid cid) const { return _gPath + Chrom::AbbrName(cid); }

	// Returns full service chrom name by chrom ID 
	const string ServName(chrid cid) const { return _sPath + Chrom::AbbrName(cid); }

	cIter begin() const { return cBegin(); }
	cIter end() const { return cEnd(); }

	//chrlen Size(const Iter::value_type& sz) const { return sz.second.Data.Real; }

	chrlen operator[] (chrid cID) const { return At(cID).Data.Real; }

	// Gets total size of genome.
	genlen GenSize() const;

#ifdef _DEBUG
	void Print() const;
#endif
};

#ifdef _ISCHIP

// 'ChromSizesExt' provides additional functionality for the ChromSizes
class ChromSizesExt : public ChromSizes
{
	mutable chrid	_treatedCnt;	// number of treated chromosomes

	// Gets chrom's effective (treated) real length: a double length for autosomes, a single somatic
	//	@it: ChromSizes iterator
	chrlen SetEffLength(cIter it) const { return Data(it).SetEffDefined(Chrom::IsAutosome(CID(it))); }

public:
	// Creates and initializes an instance.
	//	@gName: reference genome directory
	//	@customChrOpt: id of 'custom chrom' option
	//	@printMsg: true if print message about chrom.sizes generation (in case of reference genome)
	ChromSizesExt(const char* gName, BYTE customChrOpt, bool printMsg, const char* sPath)
		: _treatedCnt(0), ChromSizes(gName, customChrOpt, printMsg, sPath, true) {}

	// Gets chrom's defined effective (treated) length
	//	@it: ChromSizes iterator
	chrlen DefEffLength(cIter it) const;

	// Gets count of treated chromosomes.
	chrid TreatedCount() const { return _treatedCnt; }

	// Sets actually treated chromosomes according template and custom chrom
	//	@templ: template bed or NULL
	//	return: number of treated chromosomes
	chrid	SetTreatedChroms(bool statedAll, const Features* const templ);

	// Prints threated chroms short names
	void	PrintTreatedChroms() const;
};

#endif	// _ISCHIP

// 'ChromSeq' represented chromosome as an array of nucleotides
class ChromSeq
{
private:
	chrid	_ID;			// ID of chromosome
	char*	 _seq = NULL;	// the nucleotides buffer
	chrlen	_len;			// length of chromosome
	chrlen	_gapLen;		// total length of gaps
	Region	_effDefRgn;		// effective defined region (except 'N' at the begining and at the end)

	// Initializes instance and/or chrom's defined regions
	//	@fName: file name
	//	@rgns: chrom's defined regions: ripe or new
	//	@fill: if true fill sequence and def regions, otherwise def regions only
	//	return: true if chrom def regions are stated
	bool Init(const string& fName, ChromDefRegions& rgns, bool fill);

public:
	static bool	LetGaps;	// if true then include gaps at the edges of the ref chrom while reading
	static bool	StatGaps;	// if true count sum gaps for statistic output

	~ChromSeq() { delete[] _seq; }

	// Gets chrom legth
	chrlen Length()	const { return _len; }

	// Gets Read on position or NULL if exceeding the chrom length
	//const char* Read(chrlen pos) const { return pos > _len ? NULL : _seq + pos;	}

	//const char* Read(chrlen pos, readlen len) const { return pos + len > Length() ? NULL : _seq + pos; }

	// Gets subsequence without exceeding checking 
	const char* Seq(chrlen pos) const { return _seq + pos; }

#if defined _ISCHIP || defined _VALIGN

	// Creates a stub instance (for sampling cutting)
	//	@len: chrom length
	ChromSeq(chrlen len) : _ID(Chrom::UnID), _seq(NULL), _len(len), _gapLen(0)
	{
		_effDefRgn.Set(0, len);
	}

	// Creates and fills new instance
	ChromSeq(chrid cID, const ChromSizes& cSizes);

#endif
#ifdef _ISCHIP
	// Returns chrom ID
	chrid ID() const { return _ID; }

	// Gets count of nucleotides outside of defined region
	chrlen UndefLength() const { return Length() - _effDefRgn.Length(); }

	//float UndefLengthInPerc() const { return 100.f * (Length() - _effDefRgn.Length()) / Length(); }

	// Gets feature with defined region (without N at the begining and at the end) and score of 1
	const Featr DefRegion() const { return Featr(_effDefRgn); }

	// Gets start position of first defined nucleotide
	chrlen Start()	const { return _effDefRgn.Start; }

	// Gets end position of last defined nucleotide
	chrlen End()		const { return _effDefRgn.End; }

	// Gets total length of gaps
	chrlen GapLen()	const { return _gapLen; }

#elif defined _READDENS || defined _BIOCC

	// Creates an empty instance and fills chrom's defined regions
	//	@fName: FA file name with extension
	//	@rgns: new chrom's defined regions
	//	@minGapLen: minimal length which defines gap as a real gap
	ChromSeq(const string& fName, ChromDefRegions& rgns, short minGapLen);

#endif	// _ISCHIP

	//#if defined _FILE_WRITE && defined DEBUG 
	//	// Saves instance to file by fname
	//	void Write(const string & fname, const char *chrName) const;
	//#endif
};

#if defined _READDENS || defined _BIOCC

// 'DefRegions' represents chrom's defined regions
// initialized from ChromSizes (statically, at once)
// or from .fa files (dynamically, by request).
class DefRegions : public Chroms<Regions>
{
	ChromSizes& _cSizes;
	const chrlen	_minGapLen;	// minimal allowed length of gap
#ifdef _BIOCC
	const bool		_singleRgn = true;	// true if this instance has single Region for each chromosome
#endif


public:
	// Creates an instance by genome name, from chrom sizes file or genome.
	//	@cSizes: chrom sizes
	//	@minGapLen: minimal length which defines gap as a real gap
	DefRegions(ChromSizes& cSizes, chrlen minGapLen)
		: _cSizes(cSizes), _minGapLen(minGapLen)
	{
		Init();
	}

	void Init();

	// Returns true if regions are not initialized
	bool IsEmpty() const { return !Count(); }

	ChromSizes& ChrSizes() { return _cSizes; }

	// Gets chrom's size by chrom's iterator
	chrlen Size(cIter it)	const { return Data(it).LastEnd(); }

	// Gets chrom's size by chrom's ID
	chrlen Size(chrid cID)	const { return At(cID).Data.LastEnd(); }

	// Gets chrom regions by chrom ID; lazy for real chrom regions
	const Regions& operator[] (chrid cID);

#ifdef _BIOCC
	// Copying constructor: creates empty copy!
	DefRegions(const DefRegions& gRgns) :
		_cSizes(gRgns._cSizes), _minGapLen(gRgns._minGapLen), _singleRgn(gRgns._singleRgn) {}

	// Gets miminal size of chromosome: for represented chromosomes only
	chrlen MinSize() const;

	// Gets total genome's size: for represented chromosomes only
	genlen GenSize() const;

	// Adds chromosomes and regions without check up
	void AddChrom(chrid cID, const Regions& rgns) { AddVal(cID, rgns); }
#endif	// _BIOCC

#ifdef _DEBUG
	void Print() const;
#endif
};
#endif	// _READDENS || _BIOCC
