/**********************************************************
DataReader.cpp
Last modified: 11/22/2023
***********************************************************/

#include "DataReader.h"
#include "ChromData.h"

/************************ DataReader ************************/

// Sets the next chromosome as the current one if they are different
//	@cID: next chrom ID
//	@return: true, if new chromosome is set as current one
bool DataReader::SetNextChrom(chrid cID)
{
	if (cID == _cID || cID == Chrom::UnID)	return false;
	_cID = cID;
	return true;
}

/************************ end of DataReader ************************/

/************************ BedReader ************************/

// Reset WIG type, score index, chrom mark position offset and estimated number of lines
void BedReader::ResetWigType(FT::eType type, BYTE scoreInd, size_t cMarkPosOffset)
{
	ResetType(type);
	_scoreInd = scoreInd;
	_chrMarkPos += BYTE(cMarkPosOffset);
}

// Inserts '0' after chrom in current line and returns point to the next decl parameter if exists
//const char* BedReader::SplitLineOnChrom()
//{
//	char* line = strchr(ChromMark(), SPACE);
//	if (!line)	return NULL;
//	*line = cNULL; return line + 1;
//}

// Returns true if item contains the strand sign
//	Is invoked in the Feature constructor only.
//bool BedReader::IsItemHoldStrand() const
//{
//	if (StrandFieldInd){
//		const char s = *StrField(StrandFieldInd);
//		return s == '+' || s == '-';
//	}
//	return false;
//}

	// Checks for Fixed or Variable step wiggle type and corrects it if found
	//	@line: possible declaration line
	//	return: true if Fixed or Variable step type is specified
bool BedReader::DefineWigType(const char* line)
{
	FT::eType type = FT::eType::UNDEF;

	if (KeyStr(line, WigFixSTEP))
		ResetWigType(type = FT::eType::WIG_FIX, 0, WigFixSTEP.length() + 1);
	else if (KeyStr(line, WigVarSTEP))
		ResetWigType(type = FT::eType::WIG_VAR, 1, WigVarSTEP.length() + 1);
	if (type != FT::eType::UNDEF) {
		SetEstLineCount(type);
		RollBackRecord(TAB);					// roll back the read declaration line
		return true;
	}
	return false;
}

// Creates new instance for reading and open file; specifies the WIG type if initial type is BedGrap
//	@fName: name of file
//	@type: file type
//	@scoreNumb: number of 'score' filed (0 by default for ABED and BAM)
//	@msgFName: true if file name should be printed in exception's message
//	@abortInval: true if invalid instance should be completed by throwing exception
BedReader::BedReader(const char* fName, FT::eType type, BYTE scoreNumb, bool msgFName, bool abortInval) :
	_scoreInd(scoreNumb ? scoreNumb - 1 : 4),	// default score index for ABED and BAM
	_chrMarkPos(BYTE(strlen(Chrom::Abbr))),
	TabReader(fName, type, eAction::READ, false, msgFName, abortInval)
{
	if (type == FT::eType::ABED)
		_getStrand = [this]() { return *StrField(StrandFieldInd) == PLUS; };
	else			// for BED ignore strand to omit duplicates even when strand is defined
		_getStrand = []() { return true; };

	// ** read track definition line and clarify types
	const char* line = TabReader::GetNextLine(false);
	if (!line) ThrowExcept(Err::F_EMPTY);

	const char* line1 = KeyStr(line, "track");			// check track key
	if (line1)											// track definition line
		if (type == FT::eType::BGRAPH) {				// defined by extention
			static const char* typeBGraph = BedGraphTYPE;
			static const char* typeWiggle = WigTYPE;

			// ** define type by track definition line1
			line1 = CheckSpec(line1 + 1, "type=");
			const size_t len = strchr(line1, SPACE) - line1;	// the length of wiggle type in definition line1
			if (!len)	ThrowExcept("track type is not specified");
			if (strncmp(line1, typeBGraph, max(len, strlen(typeBGraph)))) {		// not BedGraph; be fit len == strlen(typeBGraph)
				// check for 'wiggle_0'
				if (strncmp(line1, typeWiggle, max(len, strlen(typeWiggle))))
					ThrowExcept("type '" + string(line1, len) + "' does not supported");
				// fixed or variable step
				line1 = GetNextLine(false);
				if (!DefineWigType(line1))
					ThrowExcept(string(line1) + ": absent or unknown wiggle data format");
			}
			else SetEstLineCount(type);
		}
		else SetEstLineCount();			// ordinary bed
	else								// no track definition line
		if (type == FT::eType::BGRAPH) {
			if (!DefineWigType(line)) {
				SetEstLineCount(type);	// real bedgraph
				RollBackRecord(TAB);	// roll back the read data line
			}
		}
		else SetEstLineCount();			// ordinary bed
}

/************************ end of BedReader ************************/

#ifdef _BAM
/************************ BamReader ************************/

// Creates new instance for reading and open file
//	@fName: name of file
//	@cSizes: chrom sizes to be initialized or NULL
//	@prName: true if file name should be printed in exception's message
BamReader::BamReader(const char* fName, ChromSizes* cSizes, bool prName) : _prFName(prName)
{
	_reader.Open(fName);

	// variant of estimation with max/min ~ 18
	float x = float(_reader.GetReferenceCount()) * 200000;
	_estItemCnt = ULONG(log(x * x) * 100000);
	// variant or estimation with max/min ~ 28
	//_estItemCnt = _reader.GetReferenceCount() * 200000;

#ifndef _NO_CUSTOM_CHROM
	const string headerSAM = GetHeaderText();
	if (cSizes)
		cSizes->Init(headerSAM);
	else			// validate all chroms ID via SAM header; empty functor
		Chrom::ValidateIDs(headerSAM, [](chrid, const char*) {});
#endif // _NO_CUSTOM_CHROM
}

/************************ end of BamReader ************************/
#endif	// _BAM

/************************ UniBedReader ************************/

bool UniBedReader::IsTimer = false;	// if true then manage timer by Timer::Enabled, otherwise no timer

// Returns chrom size
//	Defined in cpp because of call in template function (otherwise ''ChromSize' is no defined')
chrlen UniBedReader::ChromSize(chrid cID) const
{
	return (*_cSizes)[cID];
}

// Resets the current accounting of items
void UniBedReader::ResetChrom()
{
	_rgn0.Set();
	_issues[DUPL].Cnt += _cDuplCnt;
	_cDuplCnt = 0;
	_cCnt++;
}

// Validate item
//	@cLen: current chrom length or 0 if _cSizes is undefined
//	return: true if item is valid
bool UniBedReader::CheckItem(chrlen cLen)
{
	bool res = true;
	if (_checkSorted && _rgn.Start < _rgn0.Start)
		_file->ThrowExceptWithLineNumb("unsorted " + FT::ItemTitle(_type, false));
	if (_rgn.Invalid())
		_file->ThrowExceptWithLineNumb("'start' position is equal or more than 'end' one");
	if (cLen)				// check for not exceeding the chrom length
		if (_rgn.Start >= cLen) { _issues[STARTOUT].Cnt++; return false; }
		else if (_rgn.End > cLen) {
			_issues[ENDOUT].Cnt++;
			if (_type != FT::eType::ABED && _type != FT::eType::BAM)
				_rgn.End = cLen;
			else	return false;
		}

	_strand = _file->ItemStrand();				// the only reading strand from file
	if (_rgn0 == _rgn && _strand0 == _strand)	// duplicates
		_cDuplCnt++,
		res = _MaxDuplCnt == vUNDEF || ++_duplCnt < _MaxDuplCnt;
	else {
		_duplCnt = 0;
		_lenFreq[_file->ItemLength()]++;
		res = ChildCheckItem();					// RBed: rlen accounting; FBed: overlap check
	}
	_strand0 = _strand;
	return res;
}

// Prints count of items
//	@cnt: total count of items
//	@title: item title
void UniBedReader::PrintItemCount(ULONG cnt, const string& title)
{
	dout << SepCl;
	if (Chrom::CustomID() == Chrom::UnID)
		dout << sTotal << SPACE << cnt << SPACE << title;
	else
		dout << cnt << SPACE << title << " per " << Chrom::ShortName(Chrom::CustomID());
}

// Prints items statistics
//	@cnt: total count of items
void UniBedReader::PrintStats(ULONG cnt)
{
	PrintItemCount(cnt, FT::ItemTitle(_type, cnt != 1));
	if (cnt) {
		ULONG issCnt = 0;
		for (const Issue& iss : _issues)	issCnt += ULONG(iss.Cnt);
		if (issCnt) {
			if (_MaxDuplCnt == vUNDEF)		_issues[DUPL].Action = ACCEPT;
			else if (_MaxDuplCnt) {
				stringstream ss(" except for the first ");
				ss << _MaxDuplCnt;
				_issues[DUPL].Ext = ss.str().c_str();
			}
			_issues[OVERL].Action = GetOverlAction();
			if (_type == FT::eType::BED)	_issues[ENDOUT].Action = TRUNC;

			PrintStats(cnt, issCnt, _issues, _oinfo == eOInfo::STAT);
		}
	}
	if (!(Timer::Enabled && IsTimer))	dout << LF;
};

// Prints part number and percent of total
//	@part: part number
//	@total: total number
//	@fwidth: field width
void PrintValAndPercent(ULONG part, ULONG total, BYTE fwidth = 1)
{
	dout << setfill(SPACE) << setw(fwidth) << SPACE;
	dout << part << sPercent(part, total, 2, 0, true);
}

// Prints items statistics
//	@cnt: total count of items
//	@issCnt: count of item issues
//	@issues: issue info collection
//	@prStat: it TRUE then print issue statsistics
void UniBedReader::PrintStats(ULONG cnt, ULONG issCnt, const vector<Issue>& issues, bool prStat)
{
	static const char* sActions[] = { "accepted", "truncated", "joined", "omitted" };
	const BYTE pWidth = 4;		// padding width

	if (prStat)		dout << ", from which\n";
	for (const Issue& iss : issues)
		if (iss.Cnt) {
			if (prStat) {
				PrintValAndPercent(iss.Cnt, cnt, pWidth);
				dout << SPACE << iss.Title << SepCl << sActions[iss.Action];
				if (iss.Ext)	dout << iss.Ext;
				dout << LF;
			}
			if (iss.Action <= UniBedReader::eAction::TRUNC)
				issCnt -= iss.Cnt;
		}
	if (prStat)	dout << setw(pWidth) << SPACE << sTotal;
	else		dout << COMMA;
	dout << SPACE << sActions[0], PrintValAndPercent(cnt - issCnt, cnt);
};

UniBedReader::UniBedReader(const char* fName, const FT::eType type, ChromSizes* cSizes,
	BYTE scoreNumb, char dupLevel, eOInfo oinfo, bool prName, bool checkSorted, bool abortInval) :
	_type(type), _MaxDuplCnt(dupLevel), _checkSorted(checkSorted), _abortInv(abortInval), _oinfo(oinfo), _cSizes(cSizes)
{
	if (prName) { dout << fName; fflush(stdout); }

#ifdef _BAM
	if (type == FT::eType::BAM)
		_file = new BamReader(fName, cSizes, prName);
	else
#endif	//_BAM
		if (type <= FT::eType::ABED || type == FT::eType::BGRAPH) {
			_file = new BedReader(fName, type, scoreNumb, !prName, abortInval);
			_type = ((BedReader*)_file)->Type();	// possible change of BGRAPH with WIG_FIX or WIG_VAR
		}
		else
			Err(
#ifndef _BAM
				type == FT::eType::BAM ? "this build does not support bam files" :
#endif
				"wrong extension",
				prName ? nullptr : fName
			).Throw(abortInval);
}

UniBedReader::~UniBedReader()
{
#ifdef _BAM
	if (_type == FT::eType::BAM)
		delete (BamReader*)_file;
	else
#endif
		delete (BedReader*)_file;
}

// Returns estimated number of items
size_t UniBedReader::EstItemCount() const
{
	const size_t extCnt = _file->EstItemCount();

	//cout << LF << extCnt << TAB << _cSizes->GenSize() << TAB << Chrom::AbbrName(Chrom::CustomID()) << TAB << ChromSize(Chrom::CustomID()) <<LF;
	if (!_cSizes)	return extCnt;
	if (!Chrom::IsCustom()) {				// user stated chrom
		size_t cnt = size_t((double(extCnt) / _cSizes->GenSize()) * ChromSize(Chrom::CustomID()));
		return cnt < 2 ? extCnt : cnt;		// if data contains only one chrom, cnt can by near to 0
	}
	return extCnt;
}

/************************ end of UniBedReader ************************/

#ifdef _FEATURES
/************************ FBedReader ************************/

// Returns: true if feature is valid
bool FBedReader::ChildCheckItem()
{
	_issues[UniBedReader::eIssue::OVERL].Cnt += (_isOverlap = IsOverlap());
	return _action();
}

// Creates new instance for reading and open file
//	@fName: file name
//	@cSizes: chrom sizes
//	@scoreNmb: number of 'score' filed
//	@action: action for overlapping features
//	@prName: true if file name should be printed unconditionally
//	@abortInval: true if invalid instance should be completed by throwing exception
FBedReader::FBedReader(const char* fName, ChromSizes* cSizes,
	BYTE scoreNmb, eAction action, eOInfo oinfo, bool prName, bool abortInval) :
	_isJoin(action == eAction::JOIN),
	_overlAction(action),
	UniBedReader(fName, FT::eType::BED, cSizes, scoreNmb, 0, oinfo, prName, true, abortInval)
{
	switch (action) {
	case eAction::ACCEPT:
	case eAction::JOIN:	 _action = []() { return true; }; break;
	case eAction::OMIT:	 _action = [this]() { return !_isOverlap; }; break;
	case eAction::ABORT: _action = [this]() {
		if (_isOverlap)
			ThrowExceptWithLineNumb("overlapping features");
		return true; };
	}
}

// Returns true if features length distribution is degenerate
bool FBedReader::NarrowLenDistr() const
{
	if (!_lenFreq.size())		return false;		// no readed features
	if (_lenFreq.size() == 1)	return true;		// all features have the same length
	ULONG sum = 0;				// sum of frequencies
	for (const auto& f : _lenFreq)	sum += f.second;
	return prev(_lenFreq.cend())->second / sum > 0.9;
}

/************************ end of FBedReader ************************/
#endif	// _FEATURES

#ifndef _WIGREG

/************************ class Read ************************/

static const string TipNoFind = "Cannot find ";

#ifdef _READS
long GetNumber(const char* str, const RBedReader& file, const string& tipEnd)
{
	if (!str || !isdigit(*(++str)))
		file.ThrowExceptWithLineNumb(TipNoFind + tipEnd);
	return atol(str);
}
#endif	// _READS

readlen	Read::FixedLen;				// length of Read
const char	Read::Strands[] = { '+', '-' };

#ifdef _ISCHIP

char	Read::SeqQuality;		// the quality values for the sequence (ASCII)
bool	Read::PosInName;		// true if Read name includes a position
readlen	Read::LimitN = vUNDEF;	// maximal permitted number of 'N' in Read or vUNDEF if all
const char Read::Complements[] = { 'T',0,'G',0,0,0,'C',0,0,0,0,0,0,cN,0,0,0,0,0,'A' };
const char* Read::Title = "Read";
const char* Read::title = "read";

Read::pCopyRead Read::CopyRead[] = { &Read::Copy, &Read::CopyComplement };

//void (*Read::spCopyRead[])(const Read*, char*) = {
//	[](const Read* r, char* dst) -> void { r->Copy(dst); },
//	[](const Read* r, char* dst) -> void { r->CopyComplement(dst); }
//};

// Initializes static members
//	@len: length of Read
//	@posInName: true if Read position is included in Read name
//	@seqQual: quality values for the sequence
//	@limN: maximal permitted number of 'N'
void Read::Init(readlen len, bool posInName, char seqQual, short limN)
{
	FixedLen = len;
	PosInName = posInName;
	SeqQuality = seqQual;
	LimitN = limN;
}

// Copies complemented Read.
void Read::CopyComplement(char* dst) const
{
	const char* seq = _seq - 1;
	for (int i = Length() - 1; i >= 0; --i)
		dst[i] = Complements[(*++seq & ~0x20) - 'A'];		// any *src to uppercase
}

// Checks Read for number of 'N'
//	return:	1: NULL Read; 0: success; -1: N limit is exceeded
int Read::CheckNLimit() const
{
	if (!_seq)		return 1;
	if (LimitN != vUNDEF) {
		const char* seq = _seq - 1;
		readlen cntN = 0;
		for (int i = Length() - 1; i >= 0; --i)
			if (*++seq == cN && ++cntN > LimitN)
				return -1;
	}
	return 0;
}

// Prints Read values - parameters
//	@signOut: output marker
//	@isRVL: true if Read variable length is set
void Read::PrintParams(const char* signOut, bool isRVL)
{
	cout << signOut << Title << SepDCl;
	if (isRVL) cout << "minimum ";
	cout << "length = " << int(FixedLen);
	if (IsPosInName())	cout << SepSCl << "name includes position";
	cout << SepSCl << "N-limit" << SepCl;
	if (LimitN == readlen(vUNDEF))	cout << Options::BoolToStr(false);
	else							cout << LimitN;
	cout << LF;
}
#else

void Read::InitBase(const RBedReader& file)
{
	const Region& rgn = file.ItemRegion();
	Pos = rgn.Start;
	Len = rgn.Length();
	Strand = file.ItemStrand();
}

#ifdef _PE_READ

// PE Read constructor
Read::Read(const RBedReader& file)
{
	Numb = GetNumber(
		strrchr(file.ItemName() + 1, Read::NmNumbDelimiter),
		file,
		"number in the read's name. It should be '*.<number>'"
	);
	InitBase(file);
}

#elif defined _VALIGN

// Extended (with saved chrom & position in name) Read constructor
Read::Read(const RBedReader& file)
{
	//static const string tip1 = " in the read's name. It should be '*:<pos>.<number>'";

	const char* ss = strchr(file.ItemName() + 1, Read::NmDelimiter);
	if (ss)		ss = strstr(++ss, Chrom::Abbr);		// to be sure that 'chr' is in ss
	if (!ss)
		file.ThrowExceptWithLineNumb(TipNoFind + "chrom mark in the read's name. It should be '*chr<x>*'");
	RecCID = Chrom::ID(ss += strlen(Chrom::Abbr));
	RecPos = GetNumber(
		strchr(++ss, Read::NmPos1Delimiter),
		file,
		"position in the read's name. It should be ' * :<pos>. < number>'"
	);

	//ss = strchr(++ss, Read::NmNumbDelimiter);	// "number" position, begining with '.'
	//if (!ss || !isdigit(*(++ss)))
	//	file.ThrowExceptWithLineNumb(TipNoFind + "number" + tip1);
	//Numb = atol(ss);

	InitBase(file);
	Score = file.ItemValue();
}

void Read::Print() const
{
	dout << Pos << TAB << Strand << TAB << Chrom::AbbrName(RecCID)
		<< RecPos << TAB << setprecision(1) << Score << TAB << LF;
}

#endif	// no  _VALIGN
/************************ end of struct Read ************************/
#endif
#endif