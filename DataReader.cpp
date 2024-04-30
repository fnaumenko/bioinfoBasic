/**********************************************************
DataReader.cpp
Last modified: 04/30/2024
***********************************************************/

#include "DataReader.h"
#include "ChromData.h"

/************************ DataReader ************************/

bool DataReader::SetNextChrom(chrid cID)
{
	if (cID == _cID || cID == Chrom::UnID)	return false;
	_cID = cID;
	return true;
}

/************************ end of DataReader ************************/

/************************ BedReader ************************/

void BedReader::ResetWigType(FT::eType type, BYTE scoreInd, BYTE cMarkPosOffset)
{
	ResetType(type);
	_scoreInd = scoreInd;
	_chrMarkPos += cMarkPosOffset;
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

bool BedReader::DefineWigType(const char* line)
{
	FT::eType type = FT::UNDEF;

	if (KeyStr(line, FT::WigFixSTEP))
		ResetWigType(type = FT::WIG_FIX, 0, BYTE(FT::WigFixSTEP.length() + 1));
	else if (KeyStr(line, FT::WigVarSTEP))
		ResetWigType(type = FT::WIG_VAR, 1, BYTE(FT::WigVarSTEP.length() + 1));
	if (type != FT::UNDEF) {
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
	if (type == FT::ABED)
		_getStrand = [this]() { return *StrField(StrandFieldInd) == PLUS; };
	else		// for BED ignore strand to omit duplicates in CheckItem() even when strand is defined
		_getStrand = []() { return true; };

	// ** read track definition line and clarify types
	const char* line = TabReader::GetNextLine(false);
	if (!line) ThrowExcept(Err::F_EMPTY);

	const char* line1 = KeyStr(line, "track");			// check track key
	if (line1)											// track definition line
		if (type == FT::BGRAPH) {				// defined by extention
			// ** define type by track definition line1
			line1 = CheckSpec(line1 + 1, "type=");
			const size_t len = strchr(line1, SPACE) - line1;	// the length of wiggle type in definition line1
			if (!len)	ThrowExcept("track type is not specified");
			if (strncmp(line1, FT::BedGraphTYPE, max(len, strlen(FT::BedGraphTYPE)))) {		// not BedGraph; be fit len == strlen(typeBGraph)
				// check for 'wiggle_0'
				if (strncmp(line1, FT::WigTYPE, max(len, strlen(FT::WigTYPE))))
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
		if (type == FT::BGRAPH) {
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
	else			// validate all chroms ID via SAM header; empty function
		Chrom::ValidateIDs(headerSAM, [](chrid, const char*) {});
#endif // _NO_CUSTOM_CHROM
}

/************************ end of BamReader ************************/
#endif	// _BAM

/************************ UniBedReader ************************/

bool UniBedReader::IsTimer = false;	// if true then manage timer by Timer::Enabled, otherwise no timer

chrlen UniBedReader::ChromSize(chrid cID) const
{
	return (*_cSizes)[cID];
}

void UniBedReader::ResetChrom()
{
	_rgn0.Set();
	_issues[DUPL].Cnt += _cDuplCnt;
	_cDuplCnt = 0;
	_cCnt++;
}

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
			if (_type != FT::ABED && _type != FT::BAM)
				_rgn.End = cLen;
			else	return false;
		}

	_strand = _file->ItemStrand();				// single reading strand from file
	if (_rgn0 == _rgn && _strand0 == _strand)	// duplicates
		_cDuplCnt++,
		res = _MaxDuplCnt == BYTE_UNDEF || ++_duplCnt < _MaxDuplCnt;
	else {
		_duplCnt = 0;
		_lenFreq[_rgn.Length()]++;
		res = ChildCheckItem();					// RBed: rlen accounting; FBed: overlap check
	}
	_strand0 = _strand;
	return res;
}

void UniBedReader::PrintItemCount(size_t cnt, const string& title)
{
	dout << SepCl;
	if (cnt)
		if (Chrom::UserCID() == Chrom::UnID)
			dout << sTotal << SPACE << cnt << SPACE << title;
		else
			dout << cnt << SPACE << title << " per " << Chrom::ShortName(Chrom::UserCID());
	else
		dout << Chrom::NoChromMsg() << " in this sequence\n";
}

void UniBedReader::PrintStats(size_t cnt)
{
	PrintItemCount(cnt, FT::ItemTitle(_type, cnt != 1));
	if (cnt) {
		size_t issCnt = 0;
		for (const Issue& iss : _issues)	issCnt += iss.Cnt;
		if (issCnt) {
			if (_MaxDuplCnt == BYTE_UNDEF)	_issues[DUPL].Action = ACCEPT;
			else if (_MaxDuplCnt) {
				stringstream ss(" except for the first ");
				ss << int(_MaxDuplCnt);
				_issues[DUPL].Extra = ss.str();
			}
			_issues[OVERL].Action = GetOverlAction();
			if (_type == FT::BED)	_issues[ENDOUT].Action = TRUNC;

			PrintStats(cnt, issCnt, _issues, _oinfo == eOInfo::STAT);
		}
	}
	if (!(Timer::Enabled && IsTimer))	dout << LF;
};

// Prints part number and percent of total
//	@param part: part number
//	@param total: total number
//	@param fwidth: field width
void PrintValAndPercent(size_t part, size_t total, BYTE fwidth = 1)
{
	dout << setfill(SPACE) << setw(fwidth) << SPACE;
	dout << part << sPercent(part, total, 2, 0, true);
}

void UniBedReader::PrintStats(size_t cnt, size_t issCnt, const vector<Issue>& issues, bool prStat)
{
	static const char* sActions[] = { "accepted", "truncated", "joined", "omitted" };
	const BYTE pWidth = 4;		// padding width

	if (prStat)		dout << ", from which\n";
	for (const Issue& iss : issues)
		if (iss.Cnt) {
			if (prStat) {
				PrintValAndPercent(iss.Cnt, cnt, pWidth);
				dout << SPACE << iss.Title << SepCl << sActions[iss.Action];
				if (iss.Extra.length())	dout << iss.Extra;
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
	BYTE scoreNumb, BYTE dupLevel, eOInfo oinfo,
	bool prName, bool checkSorted, bool abortInval) :
	_type(type),
	_MaxDuplCnt(dupLevel), 
	_checkSorted(checkSorted), 
	_abortInv(abortInval), 
	_oinfo(oinfo), 
	_cSizes(cSizes)
{
	if (prName) { dout << fName; fflush(stdout); }

#ifdef _BAM
	if (type == FT::BAM)
		_file = new BamReader(fName, cSizes, prName);
	else
#endif	//_BAM
		if (type <= FT::ABED || type == FT::BGRAPH) {
			_file = new BedReader(fName, type, scoreNumb, !prName, abortInval);
			_type = ((BedReader*)_file)->Type();	// possible change of BGRAPH with WIG_FIX or WIG_VAR
		}
		else
			Err(
#ifndef _BAM
				type == FT::BAM ? "this build does not support bam files" :
#endif
				"wrong extension",
				prName ? nullptr : fName
			).Throw(abortInval);
}

UniBedReader::~UniBedReader()
{
#ifdef _BAM
	if (_type == FT::BAM)
		delete (BamReader*)_file;
	else
#endif
		delete (BedReader*)_file;
}

size_t UniBedReader::EstItemCount() const
{
	const size_t extCnt = _file->EstItemCount();

	//cout << LF << extCnt << TAB << _cSizes->GenSize() << TAB << Chrom::AbbrName(Chrom::UserCID()) << TAB << ChromSize(Chrom::UserCID()) <<LF;
	if (_cSizes && !Chrom::IsSetByUser()) {
		auto cnt = size_t((double(extCnt) / _cSizes->GenSize()) * ChromSize(Chrom::UserCID()));
		return cnt < 2 ? extCnt : cnt;		// if data contains only one chrom, cnt can by near to 0
	}
	return extCnt;
}

/************************ end of UniBedReader ************************/

//========== RBedReader

const string RBedReader::MsgNotFind = "Cannot find ";

void RBedReader::InitReadNameParser() const
{
	const char* name = ItemName();
	const char* numb = strrchr(name + 1, DOT);
	if (!numb || !isdigit(*(++numb)))
		ThrowExceptWithLineNumb(MsgNotFind + "number in the read's name. It should be '*.<number>'");
	_rNamePrefix = reclen(numb - name);
}

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
	UniBedReader(fName, FT::BED, cSizes, scoreNmb, 0, oinfo, prName, true, abortInval)
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

readlen	Read::FixedLen;				// length of Read
const char Read::Strands[] = { '+', '-' };

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

void Read::Init(readlen len, bool posInName, char seqQual, short limN)
{
	FixedLen = len;
	PosInName = posInName;
	SeqQuality = seqQual;
	LimitN = limN;
}

void Read::CopyComplement(char* dst) const
{
	const char* seq = _seq - 1;
	for (int i = Length() - 1; i >= 0; --i)
		dst[i] = Complements[(*++seq & ~0x20) - 'A'];		// any *src to uppercase
}

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

void Read::PrintParams(const char* signOut, bool isRVL)
{
	cout << signOut << Title << SepDCl;
	if (isRVL) cout << "minimum ";
	cout << "length = " << int(FixedLen);
	if (IsPosInName())	cout << SepSCl << "name includes position";
	cout << SepSCl << "N-limit" << SepCl;
	if (LimitN == readlen(vUNDEF))	cout << BoolToStr(false);
	else							cout << LimitN;
	cout << LF;
}
#else

#ifdef _PE_READ

Read::Read(const RBedReader& file) : Region(file.ItemRegion()), Strand(file.ItemStrand())
{
	if (file.IsReadNameParserUninit())
		file.InitReadNameParser();
	Numb = file.ReadNumber();
}

#elif defined _VALIGN

Read::Read(const RBedReader& file) : Region(file.ItemRegion()), Strand(file.ItemStrand()), Score(file.ItemValue())
{
	const string msgEnd = " in the read's name. It should be '*:chr<x>:<pos>.<number>'";

	if (file.IsReadNameParserUninit()) {
		const char* name = file.ItemName();
		const char* str = strchr(name + 1, Read::NmDelimiter);

		if (!str || !(str = strstr(str, Chrom::Abbr)))
			file.ThrowExceptWithLineNumb(RBedReader::MsgNotFind + "chrom mark" + msgEnd);
		file.SetReadNameParser(reclen(str + strlen(Chrom::Abbr) - name));
	}

	const char* str = file.ParsedReadName();
	RecCID = Chrom::ID(str);

	str += Chrom::MarkLength(RecCID);
	if (*str != Read::NmPos1Delimiter || !isdigit(*(++str)))
		file.ThrowExceptWithLineNumb(RBedReader::MsgNotFind + "position" + msgEnd);
	RecStart = atoui(str);
}

void Read::Print() const
{
	dout << Start << TAB << Strand << TAB << Chrom::AbbrName(RecCID)
		<< RecStart << TAB << setprecision(1) << Score << TAB << LF;
}

#endif	// no  _VALIGN
/************************ end of struct Read ************************/
#endif
#endif