/**********************************************************
ChromData.h
Provides chromosomes data functionality
Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 07/31/2024
***********************************************************/
#pragma once

#include "TxtFile.h"
#include <map>

//#define	MY_DEBUG
#define	CID(it)	(it)->first

const BYTE BYTE_UNDEF = BYTE(-1);

// 'ChromMap'
template <typename T> class ChromMap
{
public:
	typedef map<chrid, T> chrMap;
	typedef typename chrMap::iterator Iter;			// iterator
	typedef typename chrMap::const_iterator cIter;	// constant iterator

	// Returns a random-access iterator to the first element in the container
	cIter cBegin() const { return _cMap.begin(); }
	Iter Begin() { return _cMap.begin(); }

	// Returns the past-the-end iterator
	cIter cEnd() const { return _cMap.end(); }
	Iter End() { return _cMap.end(); }

	// Copies entry
	void Assign(const ChromMap& map) { _cMap = map._cMap; }

	// Returns reference to the item
	//	@param cID: chromosome's ID
	const T& At(chrid cID) const { return _cMap.at(cID); }
	T& At(chrid cID) { return _cMap.at(cID); }

	const T& operator[] (chrid cID) const { return _cMap[cID]; }
	T& operator[] (chrid cID) { return _cMap[cID]; }

	// Searches the container for a cID and returns an iterator to the element if found,
	// otherwise returns an iterator to end
	//	@param cID: chromosome's ID
	Iter GetIter(chrid cID) { return _cMap.find(cID); }
	const cIter GetIter(chrid cID) const { return _cMap.find(cID); }

	// Adds value type to the collection without checking cID
	void AddVal(chrid cID, const T& val) { _cMap[cID] = val; }
	void AddVal(chrid cID, const T&& val) { _cMap[cID] = val; }

	//// Adds value type to the collection without checking cID
	//void EmplaceVal(chrid cID, BYTE val) { _cMap.emplace(cID, T(val)); }
	////void EmplaceVal(T&&... args) { _cMap.emplace(cID, args); }

	// Removes an element from the container
	//	@param cID: chromosome's ID
	void Erase(chrid cID) { _cMap.erase(cID); }

	// Clear content
	void Clear() { _cMap.clear(); }

	// Returns true if element with cID exists in the container, and false otherwise.
	//	@param cID: chromosome's ID
	bool FindItem(chrid cID) const { return _cMap.count(cID) > 0; }

protected:
	// Returns count of elements
	size_t Count() const { return _cMap.size(); }

	// Adds class type value to the collection without checking cID. 
	// Avoids unnecessery copy constructor call
	//	@param cID: chromosome's ID
	//	@param val: class type value
	//	@returns: class type collection reference
	T& AddElem(chrid cID, const T& val) { return _cMap[cID] = val; }

	// Adds empty class type to the collection without checking cID
	//	@returns: class type collection reference
	T& AddEmptyElem(chrid cID) { return AddElem(cID, T()); }

	// Returns inner container
	const chrMap& Container() const { return _cMap; }
	chrMap& Container() { return _cMap; }

private:
	chrMap _cMap;
};

// 'ChromData' implements marked chromosome's data
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
	//	@param cit: chromosome's iterator
	bool IsTreated(typename ChromMap<ChromData<T> >::cIter cit) const { return cit->second.Treated; }

	// Returns true if chromosome by ID should be treated
	//bool IsTreated(chrid cID) const { return IsTreated(this->GetIter(cID));	}

	const T& Data(typename ChromMap<ChromData<T> >::cIter it) const { return it->second.Data; }

	T& Data(typename ChromMap<ChromData<T> >::Iter it) { return it->second.Data; }

	T& Data(chrid cID) { return this->At(cID).Data; }

	// Returns count of chromosomes
	chrid ChromCount() const { return chrid(this->Count()); }

#ifdef MY_DEBUG
	 // Gets count of treated chroms
	chrid TreatedCount() const
	{
		chrid res = 0;
		//for (const auto& c : *this) res += c.Treated;
		for (auto it = this->cBegin(); it != this->cEnd(); res += IsTreated(it++));
		return res;
	}
#endif

	// Returns true if chromosome cID exists in the container, and false otherwise.
	bool FindChrom(chrid cID) const { return this->FindItem(cID); }

	// Adds chrom size to the collection without checking cID
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
	//	@param fName: full file name
	//	@param extLen: length of file name's extention
	//	@returns: length of common prefix or -1 if there is no abbreviation chrom name in fName
	static int CommonPrefixLength(const string& fName, BYTE extLen);

	// Initializes chrom sizes from file
	void Read(const string& fName);

	// Saves chrom sizes to file
	//	@param fName: full file name
	void Write(const string& fName) const;

	// Fills external vector by chrom IDs relevant to file's names found in given directory.
	//	@param cIDs: filling vector of chrom's IDs
	//	@param gName: path to reference genome
	//	@returns: count of filled chrom's IDs
	chrid GetChromIDs(vector<chrid>& cIDs, const string& gName);

	// Initializes the paths
	//	@param gPath: reference genome directory
	//	@param sPath: service directory
	//	@param prMsg: true if print message about service fodler and chrom.sizes generation
	void SetPath(const string& gPath, const char* sPath, bool prMsg);

	// returns true if service path is defined
	bool IsServAvail() const { return _sPath.size(); }

protected:
	chrlen Length(cIter it) const { return Data(it).Real; }

public:
	const string& RefExt() const { return _ext; }

	// Creates and initializes an instance.
	//	@param gName: reference genome directory or chrom.sizes file
	//	@param prMsg: if true then print message about service fodler and chrom.sizes generation
	//	@param sPath: service directory
	//	@param checkGRef: if true then check if @gName is a ref genome dir; used in isChIP
	ChromSizes(const char* gName, bool prMsg, const char* sPath = NULL, bool checkGRef = false);

	ChromSizes() { _ext = _gPath = _sPath = strEmpty; }

	// Initializes the instance by SAM header
	void Init(const string& headerSAM);

	bool IsFilled() const { return Count(); }

	// Return true if chrom.sizes are defined explicitly, by user
	//bool IsExplicit() const { return !_gPath.length(); }		// false if path == strEmpty

	// Returns reference directory
	const string& RefPath() const { return _gPath; }

	// Returns service directory
	const string& ServPath() const { return _sPath; }

	// Returns true if the service directory coincides with the reference one
	const bool IsServAsRef() const { return _gPath == _sPath; }

	// Returns full reference chromosome's name
	//	@param cID: chromosome's ID
	const string RefName(chrid cID) const { return _gPath + Chrom::AbbrName(cID); }

	// Returns full service chromosome's name
	//	@param cID: chromosome's ID
	const string ServName(chrid cID) const { return _sPath + Chrom::AbbrName(cID); }

	cIter begin() const { return cBegin(); }
	cIter end() const { return cEnd(); }

	//chrlen Size(const Iter::value_type& sz) const { return sz.second.Data.Real; }

	chrlen operator[] (chrid cID) const { return At(cID).Data.Real; }

	// Sets all chromosomes as treated or untreated
	void TreateAll(bool treate = true);

	// Sets chromosome as treated or untreated
	void TreateChrom(chrid cID, bool treated = true) { At(cID).Treated = treated; }

	// Gets total size of genome
	genlen GenSize() const;

#ifdef MY_DEBUG
	void Print() const;
#endif
};
