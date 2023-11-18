/**********************************************************
Data.cpp (c) 2014 Fedor Naumenko (fedor.naumenko@gmail.com)
All rights reserved.
-------------------------
Last modified: 11/18/2023
-------------------------
Provides common data functionality
***********************************************************/

#include "Data.h"
#include <fstream>	// to write simple files without _FILE_WRITE

#ifdef _FEATURES
/************************ class Features ************************/

// Adds chrom to the instance
//	@cID: chrom
//	@cnt: count of chrom items
void Features::AddChrom(chrid cID, size_t cnt)
{
	if (!cnt)	return;
	const chrlen lastInd = chrlen(_items.size());
	AddVal(cID, ItemIndices(lastInd - chrlen(cnt), lastInd));
}

// treats current item
//	return: true if item is accepted
bool Features::operator()()
{
	if (_file->IsJoined()) {
		_items.back().End = _file->ItemEnd();
		return false;
	}
#ifdef _ISCHIP
	float score;
	if (_uniScore) _maxScore = score = 1;
	else {
		score = _file->ItemValue();
		if (score < 0)	_maxScore = score = _uniScore = 1;		// score is undefined in input data
		else if (score > _maxScore)	_maxScore = score;
	}
	_items.emplace_back(_file->ItemRegion(), score);
#else
	_items.emplace_back(_file->ItemRegion());
#endif
	return true;
}

// Gets the sum length of all chromosome's features
//	@cit: chromosome's iterator
chrlen Features::FeaturesLength(cIter cit) const
{
	chrlen res = 0;
	ForChrItems(cit, [&res](cItemsIter it) { res += it->Length(); });
	return res;
}

// Gets chrom's total enriched regions length:
// a double length for numeric chromosomes or a single for named.
//	@cID: chromosome's ID
//	@multiplier: 1 for numerics, 0 for nameds
//	@fLen: average fragment length on which each feature will be expanded in puprose of calculation
//	(float to minimize rounding error)
//	return: chrom's total enriched regions length, or 0 if chrom is absent
chrlen Features::EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const
{
	cIter it = GetIter(cID);
	return it != cEnd() ? EnrRegnLength(it, multiplier, fLen) : 0;
}

#ifdef _ISCHIP
// Scales defined score through all features to the part of 1.
void Features::ScaleScores()
{
	for (auto& c : Container()) {
		const auto itEnd = ItemsEnd(c.second.Data);

		for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
			it->Value /= _maxScore;		// if score is undef then it become 1
	}
}
#else	// NO _ISCHIP

// Copies feature coordinates to external DefRegions.
void Features::FillRegions(chrid cID, Regions& regn) const
{
	const auto& data = At(cID).Data;
	const auto itEnd = ItemsEnd(data);

	regn.Reserve(ItemsCount(data));
	for (auto it = ItemsBegin(data); it != itEnd; it++)
		regn.Add(*it);
}
#endif	// _ISCHIP

// Return min feature length
chrlen Features::GetMinFeatureLength() const
{
	chrlen len, minLen = CHRLEN_MAX;

	ForAllItems([&](cItemsIter it) { if ((len = it->Length()) < minLen) minLen = len; });
	return minLen;
}

// Return min distance between features boundaries
chrlen Features::GetMinDistance() const
{
	chrlen dist, minDist = CHRLEN_MAX;

	for (const auto& c : Container()) {
		const auto itEnd = ItemsEnd(c.second.Data);
		auto it = ItemsBegin(c.second.Data);

		for (chrlen end = it++->End; it != itEnd; end = it++->End)
			if ((dist = it->Start - end) < minDist)	minDist = dist;
	}
	return minDist;
}

//const chrlen UNDEFINED  = std::numeric_limits<int>::max();
#define UNDEFINED	vUNDEF

// Extends all features positions on the fixed length in both directions.
// If extended feature starts from negative, or ends after chrom length, it is fitted.
//	@extLen: distance on which start should be decreased, end should be increased
//	or inside out if it os negative
//	@cSizes: chrom sizes
//	@action: action for overlapping features
//	return: true if instance have been changed
bool Features::Extend(chrlen extLen, const ChromSizes& cSizes, UniBedInFile::eAction action)
{
	if (!extLen)	return false;
	chrlen	cRmvCnt = 0, tRmvCnt = 0;	// counters of removed items in current chrom and total removed items

	for (auto& c : Container()) {							// loop through chroms
		const chrlen cLen = cSizes.IsFilled() ? cSizes[c.first] : 0;	// chrom length
		const auto itEnd = ItemsEnd(c.second.Data);
		auto it = ItemsBegin(c.second.Data);

		it->Extend(extLen, cLen);			// first item
		cRmvCnt = 0;
		for (it++; it != itEnd; it++) {
			it->Extend(extLen, cLen);						// next item: compare to previous
			if (it->Start <= prev(it)->End)				// overlapping feature
				if (action == UniBedInFile::eAction::JOIN) {
					cRmvCnt++;
					it->Start = UNDEFINED;					// mark item as removed
					(it - cRmvCnt)->End = it->End;
				}
				else if (action == UniBedInFile::eAction::ACCEPT)
					tRmvCnt += cRmvCnt,
					cRmvCnt = 0;
				else if (action == UniBedInFile::eAction::ABORT) {
					//Err("overlapping feature with an additional extension of " + to_string(extLen)).Throw(false, true);
					dout << "overlapping feature with an additional extension of " << extLen << LF;
					return false;
				}
				else if (prev(it)->Start != UNDEFINED)		// OMIT: unmarked item
					cRmvCnt++,
					it->Start = UNDEFINED;	// mark item as removed
		}
	}
	if (cRmvCnt) {		// get rid of items marked as removed 
		vector<Featr> newItems;
		newItems.reserve(_items.size() - tRmvCnt);
		tRmvCnt = 0;
		for (auto& c : Container()) {							// loop through chroms
			ItemIndices& data = c.second.Data;
			const auto itEnd = ItemsEnd(data);
			cRmvCnt = 0;
			for (auto it = ItemsBegin(data); it != itEnd; it++)
				if (it->Start == UNDEFINED)		cRmvCnt++;	// skip removed item
				else			newItems.emplace_back(*it);
			data.FirstInd -= tRmvCnt;				// correct indexes
			data.LastInd -= (tRmvCnt += cRmvCnt);
		}
		_items.swap(newItems);
	}
	return true;
}

// Checks whether all features length exceed given length, throws exception otherwise.
//	@len: given control length
//	@lenDef: control length definition to print in exception message
//	@sender: exception sender to print in exception message
void Features::CheckFeaturesLength(chrlen len, const string& lenDef, const char* sender) const
{
	ForAllItems([&](cItemsIter it) {
		if (it->Length() < len) {
			ostringstream oss("Feature size ");
			oss << it->Length() << " is less than stated " << lenDef << SPACE << len;
			Err(oss.str(), sender).Throw();
		}
		});
}

/************************ end of class Features ************************/
#endif	// _FEATURES

/************************  class ChromSizes ************************/

// Returns length of common prefix before abbr chrom name of all file names
//	@fName: full file name
//	@extLen: length of file name's extention
//	return: length of common prefix or -1 if there is no abbreviation chrom name in fName
inline int	ChromSizes::CommonPrefixLength(const string& fName, BYTE extLen)
{
	// a short file name without extention
	return Chrom::PrefixLength(fName.substr(0, fName.length() - extLen).c_str());
}

// Initializes chrom sizes from file
void ChromSizes::Read(const string& fName)
{
	TabFile file(fName, FT::eType::CSIZE);	// file check already done

	while (file.GetNextLine()) {
		chrid cID = Chrom::ValidateIDbyAbbrName(file.StrField(0));
		if (cID != Chrom::UnID)
			AddValue(cID, ChromSize(file.LongField(1)));
	}
}

// Saves chrom sizes to file
//	@fName: full file name
void ChromSizes::Write(const string& fName) const
{
	ofstream file;

	file.open(fName.c_str(), ios_base::out);
	for (cIter it = cBegin(); it != cEnd(); it++)
		file << Chrom::AbbrName(CID(it)) << TAB << Length(it) << LF;
	file.close();
}

// Fills external vector by chrom IDs relevant to file's names found in given directory.
//	@cIDs: filling vector of chrom's IDs
//	@gName: path to reference genome
//	return: count of filled chrom's IDs
chrid ChromSizes::GetChromIDs(vector<chrid>& cIDs, const string& gName)
{
	vector<string> files;
	if (!FS::GetFiles(files, gName, _ext))		return 0;

	chrid	cid;				// chrom ID relevant to current file in files
	int		prefixLen;			// length of prefix of chrom file name
	chrid	extLen = BYTE(_ext.length());
	chrid	cnt = chrid(files.size());

	cIDs.reserve(cnt);
	sort(files.begin(), files.end());
	// remove additional names and sort listFiles
	for (chrid i = 0; i < cnt; i++) {
		if ((prefixLen = CommonPrefixLength(files[i], extLen)) < 0)		// right chrom file name
			continue;
		// filter additional names
		cid = Chrom::ValidateID(files[i].substr(prefixLen, files[i].length() - prefixLen - extLen));
		if (cid != Chrom::UnID) 		// "pure" chrom's name
			cIDs.push_back(cid);
	}
	sort(cIDs.begin(), cIDs.end());
	return chrid(cIDs.size());
}

// Initializes the paths
//	@gPath: reference genome directory
//	@sPath: service directory
//	@prMsg: true if print message about service fodler and chrom.sizes generation
void ChromSizes::SetPath(const string& gPath, const char* sPath, bool prMsg)
{
	_gPath = FS::MakePath(gPath);
	if (sPath && !FS::CheckDirExist(sPath, false))
		_sPath = FS::MakePath(sPath);
	else
		if (FS::IsDirWritable(_gPath.c_str()))
			_sPath = _gPath;
		else {
			_sPath = strEmpty;
			if (prMsg)
				Err("reference folder closed for writing and service folder is not pointed.\n").
				Warning("Service files will not be saved!");
		}
}

void  ChromSizes::SetTreatedChrom(chrid cID)
{
	if (cID == Chrom::UnID)	return;
	for (auto& c : Chroms::Container())
		c.second.Treated = c.first == cID;
}

// Creates and initializes an instance
//	@gName: reference genome directory or chrom.sizes file
//	@customChrOpt: id of 'custom chrom' option
//	@prMsg: true if print message about service fodler and chrom.sizes generation
//	@sPath: service directory
//	checkGRef: if true then check if @gName is a ref genome dir; used in isChIP
ChromSizes::ChromSizes(const char* gName, BYTE customChrOpt, bool prMsg, const char* sPath, bool checkGRef)
{
	_ext = _gPath = _sPath = strEmpty;

	Chrom::SetCustomOption(customChrOpt);
	if (gName) {
		if (FS::IsDirExist(FS::CheckedFileDirName(gName))) {	// gName is a directory
			_ext = FT::Ext(FT::eType::FA);
			SetPath(gName, sPath, prMsg);
			const string cName = _sPath + FS::LastDirName(gName) + FT::Ext(FT::eType::CSIZE);
			const bool isExist = FS::IsFileExist(cName.c_str());
			vector<chrid> cIDs;		// chrom's ID fill list

			// fill list with inizialised chrom ID and set _ext
			if (!GetChromIDs(cIDs, gName)) {				// fill list from *.fa
				_ext += ZipFileExt;			// if chrom.sizes exists, get out - we don't need a list
				if (!isExist && !GetChromIDs(cIDs, gName))	// fill list from *.fa.gz
					Err(Err::MsgNoFiles("*", FT::Ext(FT::eType::FA)), gName).Throw();
			}

			if (isExist)	Read(cName);
			else {							// generate chrom.sizes
				for (chrid cid : cIDs)
					AddValue(cid, ChromSize(FaFile(RefName(cid) + _ext).ChromLength()));
				if (IsServAvail())	Write(cName);
				if (prMsg)
					dout << FS::ShortFileName(cName) << SPACE << (IsServAvail() ? "created" : "generated") << LF,
					fflush(stdout);			// std::endl is unacceptable
			}
		}
		else {
			if (checkGRef)	Err("is not a directory", gName).Throw();
			Read(gName);		// gName is a chrom.sizes file
			_sPath = FS::DirName(gName, true);
		}
		Chrom::SetCustomID();

		SetTreatedChrom(Chrom::CustomID());
	}
	else if (sPath)
		_gPath = _sPath = FS::MakePath(sPath);	// initialized be service dir; _ext is empty!
	// else instance remains empty
}

// Initializes ChromSizes by SAM header
void ChromSizes::Init(const string& headerSAM)
{
	if (!IsFilled())
		Chrom::ValidateIDs(headerSAM,
			[this](chrid cID, const char* header) { AddValue(cID, atol(header));  }
	);
}

//void ChromSizesInit(ChromSizes* cs, const string& headerSAM)
//{
//	if (!cs->IsFilled())
//		Chrom::ValidateIDs(headerSAM,
//			[cs](chrid cID, const char* header) { cs->AddValue(cID, atol(header));  });
//};

// Gets total size of genome
genlen ChromSizes::GenSize() const
{
	if (!_gsize)
		//for(cIter it=cBegin(); it!=cEnd(); _gsize += Length(it++));
		for (const auto& c : *this)	_gsize += c.second.Data.Real;
	return _gsize;
}

#ifdef _DEBUG
void ChromSizes::Print() const
{
	//	cout << "ChromSizes: count: " << int(ChromCount()) << endl;
	//	cout << "ID\tchrom\t";
	////#ifdef _ISCHIP
	////		cout << "autosome\t";
	////#endif
	//		cout << "size\n";
	//	for(cIter it=cBegin(); it!=cEnd(); it++) {
	//		cout << int(CID(it)) << TAB << Chrom::AbbrName(CID(it)) << TAB;
	////#ifdef _ISCHIP
	////		cout << int(IsAutosome(CID(it))) << TAB;
	////#endif
	//		cout << Length(it) << TAB << LF;
	//	}

	printf("ChromSizes: count: %d\n", ChromCount());
	printf("ID  chrom   size      treated\n");
	for (const auto& c : *this)
		printf("%2d  %-8s%9d  %d\n", c.first, Chrom::AbbrName(c.first).c_str(), c.second.Data.Real, c.second.Treated);
}
#endif	// DEBUG

/************************  end of ChromSizes ************************/


#ifdef _ISCHIP
/************************  ChromSizesExt ************************/
// Gets chrom's defined effective (treated) length
//	@it: ChromSizes iterator
chrlen ChromSizesExt::DefEffLength(cIter it) const
{
	if (Data(it).Defined)	return Data(it).Defined;	// def.eff. length is initialized
	if (RefSeq::LetGaps)		return SetEffLength(it);	// initialize def.eff. length by real size
	// initialize def.eff. length by chrN.region file
	ChromDefRegions rgns(RefName(CID(it)));
	if (rgns.Empty())		return SetEffLength(it);
	return Data(it).Defined = rgns.DefLength() << int(Chrom::IsAutosome(CID(it)));
}

// Sets actually treated chromosomes according template and custom chrom
//	@templ: template bed or NULL
//	return: number of treated chromosomes
chrid ChromSizesExt::SetTreatedChroms(bool statedAll, const Features* const templ)
{
	_treatedCnt = 0;

	for (Iter it = Begin(); it != End(); it++)
		_treatedCnt +=
		(it->second.Treated = Chrom::IsCustom(CID(it))
			&& (statedAll || !templ || templ->FindChrom(CID(it))));
	return _treatedCnt;
}

// Prints threated chroms short names, starting with SPACE
void ChromSizesExt::PrintTreatedChroms() const
{
	if (TreatedCount() == ChromCount()) {
		cout << " all";
		return;
	}
	/*
	* sequential IDs printed as range: <first-inrange>'-'<last in range>
	* detached IDs or ranges are separated by comma
	*/
	chrid cID = 0, cIDlast = 0;		// current cid, last printed cid
	chrid unprintedCnt = 0;
	bool prFirst = true;	// true if first chrom in range is printed
	cIter itLast;
	auto getSep = [&unprintedCnt](chrid cnt) { return unprintedCnt > cnt ? '-' : COMMA; };
	auto printChrom = [](char sep, chrid cID) { dout << sep << Chrom::Mark(cID); };
	auto printNextRange = [&](chrid lim, chrid nextcID) {
		if (cID != cIDlast) printChrom(getSep(lim), cID);
		printChrom(COMMA, nextcID);
	};

	//== define last treated it
	for (cIter it = cBegin(); it != cEnd(); it++)
		if (IsTreated(it))	itLast = it;

	//== print treated chrom
	for (cIter it = cBegin(); it != itLast; it++) {
		if (!IsTreated(it))		continue;
		if (prFirst)
			printChrom(SPACE, cIDlast = CID(it)),
			prFirst = false;
		else
			if (CID(it) - cID > 1) {
				printNextRange(1, CID(it));
				cIDlast = CID(it);
				unprintedCnt = 0;
			}
			else
				unprintedCnt++;
		cID = CID(it);
	}

	// print last treated chrom
	cIDlast = CID(itLast);
	if (!prFirst && cIDlast - cID > 1)
		printNextRange(0, cIDlast);
	else
		printChrom(prFirst ? SPACE : getSep(0), cIDlast);
}

/************************  end of ChromSizesExt ************************/
#endif	// _ISCHIP

/************************ class RefSeq ************************/

bool RefSeq::LetGaps = true;	// if true then include gaps at the edges of the ref chrom while reading
bool RefSeq::StatGaps = false;	// if true sum gaps for statistic output

// Initializes instance and/or chrom's defined regions
//	@fName: file name
//	@rgns: chrom's defined regions: ripe or new
//	@fill: if true fill sequence and def regions, otherwise def regions only
//	return: true if chrom def regions are stated
bool RefSeq::Init(const string& fName, ChromDefRegions& rgns, bool fill)
{
	_seq = nullptr;
	bool getN = StatGaps || LetGaps || rgns.Empty();	// if true then chrom def regions should be recorded
	FaFile file(fName, rgns.Empty() ? &rgns : nullptr);

	_len = file.ChromLength();
	if (fill) {
		try { _seq = new char[_len]; }
		catch (const bad_alloc&) { Err(Err::F_MEM, fName.c_str()).Throw(); }
		const char* line = file.Line();		// First line is readed by FaFile()
		chrlen linelen;
		_len = 0;

		do	memcpy(_seq + _len, line, linelen = file.LineLength()),
			_len += linelen;
		while (line = file.NextGetLine());
	}
	else if (getN)	while (file.NextGetLine());	// just to fill chrom def regions
	file.CLoseReading();	// only makes sense if chrom def regions were filled
	_len;
	return getN;
}

#if defined _ISCHIP || defined _VALIGN

// Creates and fills new instance
RefSeq::RefSeq(chrid cID, const ChromSizes& cSizes)
{
	_ID = cID;
	ChromDefRegions rgns(cSizes.ServName(cID));	// read from file or new (empty)

	if (Init(cSizes.RefName(cID) + cSizes.RefExt(), rgns, true) && !rgns.Empty())
		_effDefRgn.Set(rgns.FirstStart(), rgns.LastEnd());
	else
		_effDefRgn.Set(0, Length());
	_gapLen = rgns.GapLen();
}

#elif defined _READDENS || defined _BIOCC

// Creates an empty instance and fills chrom's defined regions
//	@fName: FA file name with extension
//	@rgns: new chrom's defined regions
//	@minGapLen: minimal length which defines gap as a real gap
RefSeq::RefSeq(const string& fName, ChromDefRegions& rgns, short minGapLen)
{
	Init(fName, rgns, false);
	rgns.Combine(minGapLen);
}

#endif
//#if defined _FILE_WRITE && defined DEBUG
//#define FA_LINE_LEN	50	// length of wrtied lines
//
//void RefSeq::Write(const string & fName, const char *chrName) const
//{
//	FaFile file(fName, chrName);
//	chrlen i, cnt = _len / FA_LINE_LEN;
//	for(i=0; i<cnt; i++)
//		file.AddLine(_seq + i * FA_LINE_LEN, FA_LINE_LEN);
//	file.AddLine(_seq + i * FA_LINE_LEN, _len % FA_LINE_LEN);
//	file.Write();
//}
//#endif	// DEBUG

/************************ end of class RefSeq ************************/

#if defined _READDENS || defined _BIOCC

/************************ DefRegions ************************/

void DefRegions::Init()
{
	if (IsEmpty())
		if (_cSizes.IsFilled()) {
			// initialize instance from chrom sizes
			if (Chrom::IsCustom())
				for (ChromSizes::cIter it = _cSizes.cBegin(); it != _cSizes.cEnd(); it++)
					AddElem(CID(it), Regions(0, _cSizes[CID(it)]));
			else
				AddElem(Chrom::CustomID(), Regions(0, _cSizes[Chrom::CustomID()]));
			//_isEmpty = false;
		}
}

// Gets chrom regions by chrom ID; lazy for real chrom regions
const Regions& DefRegions::operator[] (chrid cID)
{
	if (FindChrom(cID))	return At(cID).Data;
	ChromDefRegions rgns(_cSizes.ServName(cID), _minGapLen);
	if (rgns.Empty())		// file with def regions doesn't exist?
	{
		//_cSizes.IsFilled();
		const string ext = _cSizes.RefExt();
		if (!ext.length())	// no .fa[.gz] file, empty service dir: _cSizes should be initialized by BAM
			return AddElem(cID, Regions(0, _cSizes[cID])).Data;
		//Err(Err::F_NONE, (_cSizes.ServName(cID) + ChromDefRegions::Ext).c_str()).Throw();
		RefSeq rs(_cSizes.RefName(cID) + ext, rgns, _minGapLen);
	}
	return AddElem(cID, rgns).Data;
}

#ifdef _BIOCC
// Gets total genome's size: for represented chromosomes only
genlen DefRegions::GenSize() const
{
	genlen gsize = 0;
	for (cIter it = cBegin(); it != cEnd(); it++)
		gsize += Size(it);
	return gsize;
}

// Gets miminal size of chromosome: for represented chromosomes only
chrlen DefRegions::MinSize() const
{
	cIter it = cBegin();
	chrlen	minsize = Size(it);
	for (it++; it != cEnd(); it++)
		if (minsize > Size(it))
			minsize = Size(it);
	return	minsize;
}
#endif	// _BIOCC

#ifdef _DEBUG
void DefRegions::Print() const
{
	cout << "DefRegions:\n";
	for (DefRegions::cIter it = cBegin(); it != cEnd(); it++)
		cout << Chrom::TitleName(CID(it))
		<< TAB << Data(it).FirstStart()
		<< TAB << Size(it) << LF;
}
#endif	// _DEBUG
/************************ DefRegions: end ************************/
#endif	// _READDENS || _BIOCC
