/**********************************************************
TxtFile.cpp
Last modified: 05/07/2024
***********************************************************/

#include "TxtFile.h"
#include <fstream>	// to write ChromDefRegions
#include <assert.h>

const BYTE TabReaderPar::BGLnLen = Chrom::MaxAbbrNameLength + 2 * 9;	// 2*pos + correction
const BYTE TabReaderPar::WvsLnLen = 9 + 3 + 2 + 25;	// pos + val + TAB + LF + correction
const BYTE TabReaderPar::WfsLnLen = 5 + 1;			// val + LF

/************************ class FT ************************/

const char* FT::bedExt = "bed";
const char* FT::wigExt = "wig";
const string FT::Interval = "interval";
const string FT::Intervals = "intervals";
const string FT::Read = "read";
const string FT::Reads = "reads";

const char* FT::BedGraphTYPE = "bedGraph";
const char* FT::WigTYPE = "wiggle_0";
const string FT::WigVarSTEP = "variableStep";
const string FT::WigFixSTEP = "fixedStep";

const FT::fTypeAttr FT::TypeAttrs[] = {
	{ "",		strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabReaderPar(1, 1) },	// undefined type
	{ bedExt,	"feature",	"features",	Mutex::eType::WR_BED,TabReaderPar(3, 6, 0, HASH, Chrom::Abbr) },	// ordinary bed
	{ bedExt,	Read,		Reads,		Mutex::eType::WR_BED,TabReaderPar(6, 6, 0, HASH, Chrom::Abbr) },	// alignment bed
	{ "sam",	strEmpty,	strEmpty,	Mutex::eType::WR_SAM,TabReaderPar(0, 0) },
	{ "bam",	Read,		Reads,		Mutex::eType::NONE,	 TabReaderPar() },
	{ wigExt,	Interval,	Intervals,	Mutex::eType::NONE,	 TabReaderPar(4, 4, TabReaderPar::BGLnLen, HASH) },	// bedgraph: Chrom::Abbr isn't specified becuase of track definition line
	{ wigExt,	Interval,	Intervals,	Mutex::eType::NONE,	 TabReaderPar(1, 1, TabReaderPar::WvsLnLen, HASH) },// wiggle_0 variable step
	{ wigExt,	Interval,	Intervals,	Mutex::eType::NONE,	 TabReaderPar(1, 1, TabReaderPar::WfsLnLen, HASH) },// wiggle_0 fixed step
	{ "fq",		Read,		Reads,		Mutex::eType::WR_FQ, TabReaderPar() },
	{ "fa",		strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabReaderPar() },
	// for the chrom.size, do not specify TabReaderPar::LineSpec to get an exception if the file is invalid
	{ "chrom.sizes",strEmpty,strEmpty,	Mutex::eType::NONE,	 TabReaderPar(2, 2, 0, cNULL) },//, Chrom::Abbr) },
	{ "region",	strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabReaderPar(2, 2) },
	{ "dist",	strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabReaderPar(1, 2) },	// required 2 data fields, but set min=1 to skip optional text line
#ifdef _ISCHIP
	{ "ini",	strEmpty,	strEmpty,	Mutex::eType::NONE,	 TabReaderPar(4, 4) },	// isChIP ini file type
#endif
};
const BYTE FT::Count = sizeof(FT::TypeAttrs) / sizeof(FT::fTypeAttr);

// Returns file type
//	@fName: file name (with case insensitive extension)
//	@isABED: if true then returns ABED in case of .bed extention
FT::eType FT::GetType(const char* fName, bool isABED)
{
	const string s_ext = FS::GetExt(fName);
	if (s_ext == "fastq")	return eType::FQ;
	const char* ext = s_ext.c_str();
	eType type = eType::UNDEF;
	for (int i = int(eType::BED); i < Count; i++)	// start from ordinary bed
		if (!_stricmp(ext, TypeAttrs[i].Extens)) {
			type = eType(i); break;
		}
	if (type == eType::BED && isABED)	type = eType::ABED;
	return type;
}

// Gets file extension, beginning at DOT and adding .gz if needed
//	@t: file type
//	@isZip: true if add ".gz"
const string FT::Ext(eType t, bool isZip)
{
	string ext = string(1, DOT) + TypeAttrs[int(t)].Extens;
	if (isZip)	ext += ZipFileExt;
	return ext;
}

/************************ end of class FT ************************/

/************************ TxtFile ************************/
const char* modes[] = { "r", "w", "a+" };
const char* bmodes[] = { "rb", "wb" };

bool TxtFile::SetBasic(const string& fName, eAction mode, void* fStream)
{
	_buff = NULL;
	_stream = NULL;
	_errCode = Err::NONE;
	_fName = fName;
	_recCnt = _currRecPos = 0;
	_buffLen = BlockSize;	// by default; can be corrected further
#ifndef _ZLIB
	if (IsZipped()) { SetError(Err::FZ_BUILD); return false; }
#endif
	// set file stream
	if (fStream)	_stream = fStream;
	else
#ifdef _ZLIB
		if (IsZipped())
			if (mode == eAction::READ_ANY)	SetError(Err::FZ_OPEN);
			else {
				if (!(_stream = gzopen(fName.c_str(), bmodes[int(mode)])))
					SetError(Err::F_OPEN);
			}
		else
#endif
			if (_stream = fopen(fName.c_str(), modes[int(mode)]))
				setvbuf((FILE*)_stream, NULL, _IONBF, 0);
			else
				SetError(Err::F_OPEN);
	return IsGood();
}

bool TxtFile::CreateIOBuff()
{
	try { _buff = new char[_buffLen]; }
	catch (const bad_alloc) { SetError(Err::F_MEM); }
	return IsGood();
}

TxtFile::TxtFile(const string& fName, eAction mode, bool msgFName, bool abortInvalid) :
	_flag(1)	// LF size set to 1
{
	string f_name(fName);	// to add '.gz' in case of zipped

	SetFlag(ZIPPED, FS::HasGzipExt(fName));
	// check for zipped option
	if (!IsZipped()
	&& mode == eAction::READ		// assumed to exist
	&& !FS::IsFileExist(f_name.c_str())) {
		f_name += ZipFileExt;
		if (!FS::IsFileExist(f_name.c_str())) {
			Err(Err::MsgNoFile(FS::ShortFileName(fName), false), FS::DirName(fName)).Throw(abortInvalid);
			return;
		}
		SetFlag(ZIPPED, true);
	}
	SetFlag(ABORTING, abortInvalid);
	SetFlag(PRNAME, msgFName);
	// try to open
	if (!SetBasic(f_name, mode, NULL))	return;
	_fSize = FS::Size(f_name.c_str());
	if (_fSize == -1)	_fSize = 0;		// new file
#ifdef _ZLIB
	else if (IsZipped()) {				// existed file
		auto size = size_t(FS::UncomressSize(f_name.c_str()));
		if (size > 0) {
			// wrong uncompressed size: increase initial size four times
			// since zip is too big to keep right size in archive
			if (size <= _fSize)	_fSize <<= 2;
			// right uncompressed size
			else				_fSize = size;
		}
		//_buffLen >>= 1;		// decrease block size twice because of allocating additional memory:
							// 2x_buffLen for writing or 3x_buffLen for reading by gzip
	}
#endif
	if (_fSize && _fSize < _buffLen)
		if (mode != eAction::WRITE)
			_buffLen = bufflen(_fSize + 1);	// for reading
		else if (_fSize * 2 < _buffLen)
			_buffLen = bufflen(_fSize * 2);	// for writing: increase small buffer for any case

	if (!CreateIOBuff())	return;

#ifdef ZLIB_NEW
	if (IsZipped() && gzbuffer((gzFile)_stream, _buffLen) == -1)
		SetError(Err::FZ_MEM);
#endif
}

#ifdef _MULTITHREAD
TxtFile::TxtFile(const TxtFile& file) : _flag(file._flag)
{
	if (!SetBasic("Clone " + file._fName, eAction::WRITE, file._stream))	return;
	RaiseFlag(CLONE);
	RaiseFlag(MTHREAD);
	file.RaiseFlag(MTHREAD);
	CreateIOBuff();
}
#endif

TxtFile::~TxtFile()
{
	if (IsClone())	return;
	if (_buff)		delete[] _buff;
	if (_stream &&
#ifdef _ZLIB
		IsZipped() ? gzclose((gzFile)_stream) :
#endif
		fclose((FILE*)_stream))
		SetError(Err::F_CLOSE);
	//cout << "Free " << _fName << "\trecords = " << _recCnt << endl;
}

void TxtFile::SetError(Err::eCode code, const string& senderSpec, const string& spec) const
{
	_errCode = code;
	if (IsFlag(ABORTING))
		Err(code, (CondFileName() + senderSpec).c_str(), spec.length() ? spec.c_str() : NULL).Throw();
}

/************************ TxtFile: end ************************/

/************************ TxtReader ************************/

TxtReader::TxtReader(const string& fName, eAction mode,
	BYTE cntRecLines, bool msgFName, bool abortInvalid) :
	_recLineCnt(cntRecLines),
	TxtFile(fName, mode, msgFName, abortInvalid)
{
	if (Length() && ReadBlock(0) >= 0)			// read first block
		_linesLen = new reclen[cntRecLines];	// nonempty file: set lines buffer
	else
		RaiseFlag(ENDREAD);						// empty file
}

int TxtReader::ReadBlock(const bufflen offset)
{
	bufflen readLen;
#ifdef _ZLIB
	if (IsZipped()) {
		int len = gzread((gzFile)_stream, _buff + offset, _buffLen - offset);
		if (len < 0) { SetError(Err::F_READ); return -1; }
		readLen = bufflen(len);
	}
	else
#endif
	{
		readLen = bufflen(fread(_buff + offset, sizeof(char), _buffLen - offset, (FILE*)_stream));
		if (readLen != _buffLen - offset && !feof((FILE*)_stream))
		{ SetError(Err::F_READ); return -1; }
	}
	_readedLen = readLen + offset;
	//#ifdef ZLIB_OLD
		//if( _readTotal + _readedLen > _fSize )
		//	_readedLen = _fSize - _readTotal;
		//_readTotal += _readedLen;
	//#endif
	if (_readedLen != _buffLen && !IsZipped() && ferror((FILE*)_stream))
	{ SetError(Err::F_READ); return -1; }
	_currRecPos = 0;
	return int(_readedLen);
}

bool TxtReader::CompleteBlock(bufflen currLinePos, bufflen blankLineCnt)
{
	if (_readedLen != _buffLen)		// final block, normal completion
	{ RaiseFlag(ENDREAD); return true; }
	if (!currLinePos)				// block is totally unreaded
	{ SetError(Err::F_BIGLINE); RaiseFlag(ENDREAD); return true; }

	blankLineCnt = _readedLen - _currRecPos - blankLineCnt;	// now length of unreaded remainder
	// move remainder to the beginning of buffer; if length of unreaded remain = 0, skip moving
	memmove(_buff, _buff + _currRecPos, blankLineCnt);
	_recLen = 0;
	if (!ReadBlock(blankLineCnt))
	{ RaiseFlag(ENDREAD); return true; };

	return false;
}

// Reads record
//	return: pointer to line or NULL if no more lines
const char* TxtReader::GetNextRecord()
{
	// Sets _currLinePos to the beginning of next non-empty line inread/write buffer
	// GetNextRecord(), GetRecordN() and GetRecordTab() are quite similar,
	// but they are separated for effectiveness
	if (IsFlag(ENDREAD))	return NULL;

	bufflen i, blanklCnt = 0;			// counter of empty lines
	bufflen currPos = _currRecPos;		// start position of current readed record
	char* buf;

	_recLen = 0;
	for (BYTE rec = 0; rec < _recLineCnt; rec++)
		for (buf = _buff + (i = currPos);; buf++, i++) {
			if (*buf == LF) {
				if (i == currPos) {				// LF marker is first in line
					++currPos; ++blanklCnt;		// skip empty line
					continue;
				}
				if (IsLFundef())	SetLF(*(buf - 1));		// define LF size
			lf:				_recLen += (_linesLen[rec] = ++i - currPos);
				currPos = i;
				break;							// mext line in a record
			}
			if (i >= _readedLen) {	// check for oversize buffer
				//cout << ">>> " << i << TAB << currPos << TAB << _readedLen << TAB << _buffLen << LF;
				if (_readedLen != _buffLen && i > currPos)	// last record does not end with LF
					goto lf;
				if (CompleteBlock(currPos, blanklCnt))	return NULL;
				currPos = blanklCnt = i = rec = 0;
				buf = _buff;
			}
		}
	_currRecPos = currPos;			// next record position
	_recCnt++;
	return RealRecord();
}

// Reads N-controlled record
//	@counterN: counter of 'N'
//	return: pointer to line or NULL if no more lines
const char* TxtReader::GetNextRecord(chrlen& counterN)
{
	if (IsFlag(ENDREAD))	return NULL;

	bufflen i, blanklCnt = 0,			// counter of empty lines
			currPos = _currRecPos;		// start position of current readed record
	chrlen cntN = 0;				// local counter of 'N'
	char* buf;

	_recLen = 0;
	for (BYTE rec = 0; rec < _recLineCnt; rec++)
		for (buf = _buff + (i = currPos);; buf++, i++) {
			if (*buf == cN)		cntN++;
			else if (*buf == LF) {
				if (i == currPos) {				// LF marker is first in line
					currPos++; blanklCnt++;		// skip empty line
					continue;
				}
				if (IsLFundef())	SetLF(*(buf - 1));		// define LF size
			lf:				_recLen += (_linesLen[rec] = ++i - currPos);
				currPos = i;
				break;
			}
			if (i >= _readedLen) {	// check for oversize buffer
				if (_readedLen != _buffLen && i > currPos)	// last record does not end with LF
					goto lf;
				if (CompleteBlock(currPos, blanklCnt))	return NULL;
				currPos = blanklCnt = i = cntN = rec = 0;
				buf = _buff;
			}
		}
	_currRecPos = currPos;			// next record position
	_recCnt++;
	counterN += cntN;
	return RealRecord();
}

// Reads tab-controlled record
//	@tabPos: TAB's positions array that should be filled
//	@cntTabs: maximum number of TABS in TAB's positions array
//	return: pointer to line or NULL if no more lines
char* TxtReader::GetNextRecord(short* const tabPos, const BYTE tabCnt)
{
	if (IsFlag(ENDREAD))	return NULL;

	bufflen i, blanklCnt = 0,			// counter of empty lines
			currPos = _currRecPos;		// start position of current readed record
	BYTE tabInd = 1;				// index of TAB position in tabPos
	char* buf;

	_recLen = 0;
	for (BYTE rec = 0; rec < _recLineCnt; rec++)
		for (buf = _buff + (i = currPos);; buf++, i++) {
			if (*buf == TAB) {
				if (tabInd < tabCnt)
					tabPos[tabInd++] = short(i + 1 - currPos);
			}
			else if (*buf == LF) {
				if (i == currPos) {				// LF marker is first in line
					currPos++; blanklCnt++;		// skip empty line
					continue;
				}
				if (IsLFundef())	SetLF(*(buf - 1));		// define LF size
			lf:				_recLen += (_linesLen[rec] = ++i - currPos);
				currPos = i;
				break;
			}
			if (i >= _readedLen) {	// check for oversize block
				if (_readedLen != _buffLen && i > currPos)	// last record does not end with LF
					goto lf;
				if (CompleteBlock(currPos, blanklCnt))	return NULL;
				currPos = blanklCnt = i = rec = 0;
				buf = _buff;
				tabInd = 1;
			}
		}
	_currRecPos = currPos;			// next record position
	_recCnt++;
	return RealRecord();
}

void TxtReader::RollBackRecord(char sep)
{
	for (auto i = _recLen; i > 1; --i)
		if (!_buff[_currRecPos - i])
			_buff[_currRecPos - i] = sep;	// return back TAB in the last record instead of inserted '0'
	_buff[_currRecPos - 1] = LF;			// return back LineFeed in the end of last record instead of inserted '0'
	_currRecPos -= _recLen;
	_recLen = 0;
	_recCnt--;
}

const string TxtReader::LineNumbToStr(Err::eCode code, BYTE lineInd) const
{
	_errCode = code;
	ostringstream ss;
	if (IsFlag(PRNAME))	ss << FileName();
	ss << ": line " << LineNumber(lineInd);
	return ss.str();
}

/************************ TxtReader: end ************************/

#ifdef _TXT_WRITER

/************************ TxtWriter ************************/

bool TxtWriter::Zipped;			// true if file should be zippped

//TxtWriter::fAddChar TxtWriter::fLineAddChar[] = {	// 'Add delimiter' methods
//	&TxtWriter::AddCharEmpty,	// empty method
//	&TxtWriter::AddDelim		// adds delimiter and increases current position
//};

bool TxtWriter::CreateLineBuff(reclen len)
{
	try { _lineBuff = new char[_lineBuffLen = len]; }
	catch (const bad_alloc) { SetError(Err::F_MEM); };
	return IsGood();
}

void TxtWriter::EndRecordToIOBuff(bufflen len)
{
	_currRecPos += len;
	_buff[_currRecPos++] = LF;
#ifdef _MULTITHREAD
	if (IsFlag(MTHREAD) && !IsClone())
		InterlockedIncrement(_totalRecCnt);		// atomic increment _recCnt
	else
#endif
		_recCnt++;
}

#ifdef _MULTITHREAD
TxtWriter::TxtWriter(const TxtWriter& file) :
	_lineBuffOffset(file._lineBuffOffset),
	_delim(file._delim),
	TxtFile(file)
{
	if (!CreateLineBuff(file._lineBuffLen))	return;
	memcpy(_lineBuff, file._lineBuff, _lineBuffLen);
	_totalRecCnt = &file._recCnt;
}
#endif

TxtWriter::~TxtWriter()
{
	if (_currRecPos)	Write();
	delete[] _lineBuff;
}

void TxtWriter::LineAddChar(char ch, bool addDelim)
{
	LineAddChar(ch);
	LineAddDelim(addDelim);
}

reclen TxtWriter::LineAddChars(const char* src, reclen len, bool addDelim)
{
	memcpy(_lineBuff + _lineBuffOffset, src, len);
	_lineBuffOffset += len;
	LineAddDelim(addDelim);
	return _lineBuffOffset;
}

void TxtWriter::LineAddInt(LLONG v, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%lld", v);
	LineAddDelim(addDelim);
}

void TxtWriter::LineAddInts(ULONG v1, ULONG v2, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%u%c%u", v1, _delim, v2);
	LineAddDelim(addDelim);
}

void TxtWriter::LineAddUInts(chrlen v1, chrlen v2, chrlen v3, bool addDelim)
{
	_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, "%u%c%u%c%u", v1, _delim, v2, _delim, v3);
	LineAddDelim(addDelim);
}

// Adds floating point value to the current position of the line write buffer,
//	adds delimiter after value and increases current position.
//	@val: value to be set
//	@ndigit: number of digits to generate
//	@addDelim: if true then adds delimiter and increases current position
//void TxtWriter::LineAddFloat(float val, BYTE ndigit, bool addDelim)
//{
//	// double val, because casting int to float is not precise, f.e. (float)61342430 = 61342432
//	_gcvt(val, ndigit, _lineBuff + _lineBuffOffset);
//	_lineBuffOffset += reclen(strlen(_lineBuff + _lineBuffOffset));
//	LineAddDelim(addDelim);
//}

void TxtWriter::SetFloatFractDigits(BYTE digitsCnt)
{
	assert(digitsCnt < 10);
	_formatFloat[0][3] = std::to_string(digitsCnt)[0];	// "%4.2f"
	_floatBuffLen = 4 + digitsCnt;
}

void TxtWriter::LineAddFloat(float val, bool addDelim)
{
	_lineBuffOffset += std::snprintf(_lineBuff + _lineBuffOffset, _floatBuffLen, _formatFloat[val == 0].c_str(), val);
	LineAddDelim(addDelim);
}

void TxtWriter::LineAddSingleFloat(float val)
{
	LineAddFloat(val, false);
	LineToIOBuff();
}

void TxtWriter::LineToIOBuff(reclen offset)
{
	RecordToIOBuff(_lineBuff, _lineBuffOffset);
	_lineBuffOffset = offset;
}

void TxtWriter::RecordToIOBuff(const char* src, bufflen len)
{
	if (_currRecPos + len + 1 > _buffLen)	// write buffer to file if it's full
		Write();
	memcpy(_buff + _currRecPos, src, len);
	EndRecordToIOBuff(len);
}

void TxtWriter::StrToIOBuff(const string&& str)
{
	memmove(_buff + _currRecPos, str.c_str(), str.length());
	EndRecordToIOBuff(bufflen(str.length()));
}

void TxtWriter::SetLineBuff(reclen len)
{
	CreateLineBuff(len);
	memset(_lineBuff, _delim, len);
}

// Adds to line 2 int values and adds line to the IO buff.
//void TxtWriter::WriteLine(ULONG val1, ULONG val2)
//{
//	LineAddInts(val1, val2, false);
//	LineToIOBuff();
//}

// Adds to line C-string with delimiter, and adds line to the IO buff.
//void TxtWriter::WriteLine(const char* str)
//{
//	LineAddStr(str, strlen(str));
//	LineToIOBuff();
//}

// Adds to line string with delimiter and int value, and adds line to the IO buff.
//void TxtWriter::WriteLine(const string& str, int val)
//{
//	LineAddStr(str);
//	LineAddInt(val, false);
//	LineToIOBuff();
//}

void TxtWriter::Write() const
{
#ifdef _MULTITHREAD
	const bool lock = IsFlag(MTHREAD) && Mutex::IsReal(_mtype);
	if (lock)	Mutex::Lock(_mtype);
#endif
	size_t res =
#ifdef _ZLIB
		IsZipped() ?
		gzwrite((gzFile)_stream, _buff, (unsigned int)_currRecPos) :
#endif
		fwrite(
			_buff,
			1,
			_currRecPos,
			(FILE*)_stream);
	if (res == _currRecPos)	_currRecPos = 0;
	else { SetError(Err::F_WRITE); /*cout << "ERROR!\n";*/ }
#ifdef _MULTITHREAD
	if (lock) {
		if (IsClone())
			//InterlockedExchangeAdd(_totalRecCnt, _recCnt);
			*_totalRecCnt += _recCnt,
			_recCnt = 0;
		Mutex::Unlock(_mtype);
	}
#endif
}

//bool	TxtWriter::AddFile(const string fName)
//{
//	TxtFile file(fName, *this);
//	while( file.ReadBlock(0) ) {
//		_currRecPos = file._readedLen;
//		if( !Write() )	return false;
//	}
//	return true;
//}

/************************ class TxtWriter: end ************************/

#endif	// _TXT_WRITER

#ifndef _FQSTATN

/************************ class TabReader ************************/

bool TabReader::IsFieldValid(BYTE ind) const
{
	if (_fieldPos[ind] == vUNDEF) {		// vUNDEF was set in buff if there was no such field in the line
	//|| !StrField(ind)[0]) {			// empty field
		if (ind < FT::FileParams(_fType).MinFieldCnt)
			Err(Err::TF_FIELD, LineNumbToStr(Err::TF_FIELD).c_str()).Throw();
		return false;
	}
	return true;
}

void TabReader::Init(FT::eType type, bool isEstLineCnt)
{
	if (!_fieldPos)
		_fieldPos = new short[FT::FileParams(type).MaxFieldCnt];
	_lineSpecLen = FT::FileParams(type).LineSpecLen();

	// estimate line cnt
	if (isEstLineCnt)
		if (FT::FileParams(type).AvrLineLen)
			SetEstLineCount(type);
		else if (GetNextLine(true))
			SetEstLineCount();
		else _estLineCnt = 0;
	//cout << "  file estLen: " << _estLineCnt << LF;
}

void TabReader::SetEstLineCount()
{
	_estLineCnt = Length() / RecordLength();
	RollBackLine();
}

void TabReader::ResetType(FT::eType type)
{
	if (FT::FileParams(_fType).MaxFieldCnt < FT::FileParams(type).MaxFieldCnt)
		Release();
	Init(_fType = type, true);
}

const char* TabReader::KeyStr(const char* str, const string& key)
{
	if (str) {
		const char* strKey = strstr(str, key.c_str());
		if (strKey)	return strKey + key.length();
	}
	return NULL;
}

const char* TabReader::CheckSpec(const char* str, const string& key)
{
	const char* strKey = KeyStr(str, key);
	if (!strKey)
		ThrowExcept(LineNumbToStr() + ": absent or wrong '" + key + "' key");
	return strKey;
}

// Skip commented lines and returns estimated number of uncommented lines
//ULONG TabReader::GetUncommLineCount()
//{
//	const char*	line;
//	const char comment = FT::FileParams(_fType).Comment;
//	ULONG	cnt = 0;
//	
//	for(USHORT pos=0; line = GetNextRecord(); pos=0) {
//		while( *(line+pos)==SPACE )	pos++;		// skip blanks at the beginning of line
//		if (*(line+pos) != comment)	break;		// skip comment line
//	}
//	if(line) {
//		cnt = ULONG(Length() / RecordLength());
//		RollBackLastRecord();
//	}
//	return cnt;
//}

// Reads first line and set it as current.
// Throws exception if invalid
//	@cntLines: returned estimated count of lines.
//	It works properly only if lines are sorted by ascending, f.e. in sorted bed-files.
//	return: current line
//const char* TabReader::GetFirstLine(ULONG& cntLines)
//{
//	cntLines = 0;
//	if(GetNextLine() )
//		cntLines = (ULONG)(Length() / RecordLength());	
//	return _currLine;
//}

const char* TabReader::GetNextLine(bool checkTab)
{
	const TabReaderPar par = FT::FileParams(_fType);

	// fill _fieldPos 0 to check if all fields will be initialize
	// _fieldPos keeps start positions (next after TAB) for each field in line.
	// First field's position usually is 0, but if there are blanks in the beginning of the line,
	// it will be in first unblank position.
	// GetNextRecord(..) fills _fieldPos beginning from second field.
	memset(_fieldPos, vUNDEF, sizeof(short) * par.MaxFieldCnt);

	char* line = GetNextRecord(_fieldPos, par.MaxFieldCnt);	// a bit faster than using _currLine
	if (line) {
		USHORT	currPos = 0;	// a bit faster than using _currPos in heep

		// skip blanks at the beginning of line
		while (line[currPos] == SPACE)	currPos++;
		// skip comment line or line without specifier
		if (*(line + currPos) == par.Comment
			|| (par.LineSpec && memcmp(line + currPos, par.LineSpec, _lineSpecLen)))
			return GetNextLine(checkTab);

		_fieldPos[0] = currPos;		// set start position of first field
#ifdef _WIG_READER
		if (checkTab || isdigit(line[0]))
#endif // _WIG_READER
			for (BYTE i = 1; i < par.MaxFieldCnt; line[_fieldPos[i++] - 1] = cNULL)	// replace TABs by 0
				if (_fieldPos[i] == vUNDEF)		// numbers of TABs is less than number of fields
					if (i >= par.MinFieldCnt)				// less than number of optional fields
						break;
					else {									// less than number of required fields
						SetError(Err::TF_FIELD, LineNumbToStr(),
							": " + to_string(i) + " against " + to_string(par.MinFieldCnt) + "; wrong format?");
						return _currLine = NULL;
					}
		line[RecordLength() - 1] = cNULL;	// replace '\n' by 0

		//cout << this->FileName() << "\trecLen = " << RecordLength() << LF;

		// this version when _fieldPos is not initialized in TxtFile::GetNextRecord(): a bit more slower
		//for(BYTE i=1; i<_cntFields; i++) {
		//	// search for the next TAB
		//	for(; currPos<_recLen; currPos++)
		//		if(line[currPos]==TAB)
		//			break;
		//	line[currPos] = cNULL;		// replace TABs by 0
		//	_fieldPos[i] = ++currPos;
		//}
	}
	return _currLine = line;
}

void TabReader::InitRegion(BYTE fInd, Region& rgn) const
{
	const char* str = StrField(fInd);
	rgn.Start = atoui_by_ref(str);
	rgn.End = atoui_by_ref(++str);
}


/************************ end of TabReader ************************/

/************************ ChromDefRegions ************************/

const string ChromDefRegions::Ext = ".region";	// regions file extension

bool ChromDefRegions::Combiner::ExceptRegion(Region& rgn)
{
	if (_rgn.Empty()) { _rgn = rgn; return false; }	// first income region
	if (rgn.Start - _rgn.End < _gapLen) { _rgn.End = rgn.End; return false; }
	Region	next = rgn;
	rgn = _rgn;		// returned region
	_rgn = next;
	// swap ???
	return true;
}

ChromDefRegions::ChromDefRegions(const string& fName, chrlen minGapLen) : _gapLen(0), _new(true)
{
	_fName = fName + Ext;
	if (FS::IsFileExist(_fName.c_str())) {
		TabReader file(_fName, FT::RGN);
		Combiner comb(minGapLen);
		Region rgn;

		Reserve(file.EstLineCount());
		if (file.GetNextLine()) {
			_gapLen = file.UIntField(1);	// first line contains total length of gaps
			while (file.GetNextLine()) {
				file.InitRegion(0, rgn);
				if (!minGapLen || comb.ExceptRegion(rgn))
					Add(rgn);
			}
			if (minGapLen)	Add(comb.LastRegion());
			_new = false;
		}
	}
	else	Reserve(DefCapacuty);
}

void ChromDefRegions::AddRegion(const Region& rgn, chrlen minGapLen)
{
	if (rgn.End && minGapLen && rgn.Length())	// current gap is closed
		if (_regions.size() && rgn.Start - _regions.back().End < minGapLen)
			_regions.back().End = rgn.End;		// pass minimal allowed gap
		else
			_regions.push_back(rgn);			// add new def-region
}

void ChromDefRegions::Write() const
{
	if (!_new || FS::IsShortFileName(_fName))	return;
	ofstream file;

	file.open(_fName.c_str(), ios_base::out);
	file << "SumGapLen:" << TAB << _gapLen << LF;
	for (Iter it = Begin(); it != End(); it++)
		file << it->Start << TAB << it->End << LF;
	file.close();
	_new = false;
}

#if defined _READDENS || defined _BIOCC
void ChromDefRegions::Combine(chrlen minGapLen)
{
	Region rgn;
	Regions rgns;	// new combined regions
	Combiner comb(minGapLen);

	rgns.Reserve(Count());
	for (Iter it = Begin(); it != End(); it++) {
		rgn.Set(it->Start, it->End);
		if (comb.ExceptRegion(rgn))	rgns.Add(rgn);
	}
	rgns.Add(comb.LastRegion());
	Copy(rgns);
}
#endif

/************************ end of ChromDefRegions ************************/

/************************ class FaReader ************************/

void FaReader::DefRgnMaker::AddGap(chrlen start, chrlen len)
{
	_defRgn.End = (start += _currPos);
	_defRgns.AddRegion(_defRgn, _minGapLen);
	_defRgn.Start = start + len;
	_defRgns.IncrGapLen(len);
}

void FaReader::DefRgnMaker::CloseAddGaps(chrlen cLen)
{
	_defRgn.End = cLen;
	_defRgns.AddRegion(_defRgn, _minGapLen);
	_defRgns.Write();
}

void FaReader::CountN(chrlen startPos, reclen NCnt)
{
	if (NCnt == 1)	return;

	const char* line = Line();
	const reclen lineLen = LineLength();
	reclen nCnt = 0, iN;				// local count N, first N index
	bool firstN = true;

	for (chrlen i = startPos; i < lineLen; i++)
		if (line[i] == cN) {					// some 'N' subsequence begins?
			nCnt++;
			if (firstN) { iN = i; firstN = false; }
		}
		else if (nCnt)						// is current subsequence finished?
			if (nCnt == 1) {					// was it a single 'N'?
				nCnt = 0; NCnt--; firstN = true;
			}
			else {							// continuous 'N' sequence?
				_rgnMaker->AddGap(iN, nCnt);	// add complete 'N' subsequence
				CountN(i++, NCnt - nCnt);	// search in the rest of line
				break;
			}
}

const char* FaReader::GetLineWitnNControl()
{
	chrlen nCnt = 0;	// number of 'N' in line
	const char* line = GetNextRecord(nCnt);
	if (line) {
		const chrlen len = LineLength();
		if (nCnt)
			if (nCnt == len)	_rgnMaker->AddGap(0, nCnt);	// the whole line is filled by 'N'
			else			CountN(0, nCnt);
		_rgnMaker->AddLineLen(len);
	}
	return line;
}

FaReader::FaReader(const string& fName, ChromDefRegions* rgns) : TxtReader(fName, eAction::READ, 1)
{
	if (rgns) {
		_rgnMaker.reset(new DefRgnMaker(*rgns, 2));
		_pGetLine = &FaReader::GetLineWitnNControl;
	}
	else _pGetLine = &FaReader::GetNextRecord;

	chrlen len = chrlen(Length());
	if (NextGetLine()[0] == FaComment) {		// is first line a header?
		len -= chrlen(RecordLength());
		if (_rgnMaker)	_rgnMaker->RemoveLineLen(LineLength());
		NextGetLine();
	}
	// set chrom length
	_cLen = len - chrlen(LFSize() *
		(len / chrlen(RecordLength()) +		// amount of LF markers for whole lines
			bool(len % chrlen(RecordLength()))));	// LF for part line
}

/************************ end of class FaReader ************************/

#endif	// no _FQSTATN

// Creates new instance with read buffer belonges to aggregated file: constructor for concatenating.
//	For reading only.
//	@fName: valid full name of file
//	@file: file that is aggregated
//TxtFile::TxtFile(const string & fName, const TxtFile& file) :
//	_flag(file._flag),
//	_recLineCnt(file._recLineCnt),
//	_buffLineLen(0)
//{
//	if( !SetBasic(fName, READ, NULL) )	return;
//	_buffLen = file._buffLen;
//	_buff = file._buff;
//	RaiseFlag(CONSTIT);	
//}