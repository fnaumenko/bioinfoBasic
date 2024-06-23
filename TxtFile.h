/**********************************************************
TxtFile.h
Provides read|write basic bioinfo text files functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 06/23/2024
***********************************************************/
#pragma once

#include "common.h"

// Number of basics file's reading|writing buffer blocks.
// Should be less than 2047 because of ULONG type of block size variable.
// Otherwise the behaviour is unpredictable.

typedef uint16_t reclen;	// type: length of record

// 'TabReaderPar' keeps basic parameters for TabReader
struct TabReaderPar
{
	static const BYTE BGLnLen;	// predefined BedGraph line length in purpose of reserving container size
	static const BYTE WvsLnLen;	// predefined WIG variable step line length in purpose of reserving container size
	static const BYTE WfsLnLen;	// predefined WIG fixed step line length in purpose of reserving container size

	const BYTE	MinFieldCnt;	// minimum number of feilds in file line; checked during initialization
	const BYTE	MaxFieldCnt;	// maximum possible number of feilds in data line; checked by a call
	const BYTE	AvrLineLen;		// average line length; to reserve container size
	const char	Comment;		// char indicates that line is comment
	const char* LineSpec;		// substring on which each data line is beginning

	TabReaderPar() : MinFieldCnt(0), MaxFieldCnt(0), AvrLineLen(0), Comment(cNULL), LineSpec(NULL) {}

	TabReaderPar(BYTE minTabCnt, BYTE maxTabCnt, BYTE avrLineLen = 0, char comm = HASH, const char* lSpec = NULL) :
		MinFieldCnt(minTabCnt),
		MaxFieldCnt(maxTabCnt < minTabCnt ? minTabCnt : maxTabCnt),
		AvrLineLen(avrLineLen),
		Comment(comm),
		LineSpec(lSpec)
	{}

	reclen LineSpecLen() const { return reclen(LineSpec ? strlen(LineSpec) : 0); }
};

// 'File Type' implements bioinformatics file type routines 
static class FT
{
	struct fTypeAttr {
		const char*	 Extens;		// file extension
		const string Item;			// item title
		const string ItemPl;		// item title in plural
		Mutex::eType MtxType;		// type of Mutex lock
		TabReaderPar FileParam;		// TabReader parameters, defined feilds
	};
	static const char*	bedExt;
	static const char*	wigExt;
	static const string Interval;
	static const string Intervals;
	static const string Read;
	static const string Reads;
	static const fTypeAttr	TypeAttrs[];
	static const BYTE	Count;

public:
	// bioinfo file types
	enum /*class*/ eType {
		UNDEF,	// undefined type
		BED,	// ordinary bed
		ABED,	// alignment bed
		SAM,	// sam
		BAM,	// bam
		BGRAPH,	// bedgraph
		WIG_VAR,// wiggle_0 variable step
		WIG_FIX,// wiggle_0 fixed step
		FQ,		// fastQ
		FA,		// fasta
		CSIZE,	// chrom sizes
		RGN,	// my own region type
		DIST,	// my own distribution type
#ifdef _ISCHIP
		INI		// isChIP ini file type
#endif
	};

	static const char* BedGraphTYPE;
	static const char* WigTYPE;
	static const string WigVarSTEP;
	static const string WigFixSTEP;

	// Gets file type
	//	@fName: file name (with case insensitive extension)
	//	@isABED: if true then returns ABED in case of .bed extention
	static eType GetType(const char* fName, bool isABED = false);

	// Gets file extension, beginning at DOT and adding .gz if needed
	//	@t: file type
	//	@isZip: true if add ".gz"
	static const string Ext(eType t, bool isZip = false);

	// Gets an item's title
	//	@t: file type
	//	@pl: true if plural form
	static const string& ItemTitle(eType t, bool pl = true)
	{
		return pl ? TypeAttrs[int(t)].ItemPl : TypeAttrs[int(t)].Item;
	}

	// Gets TabReader params
	//	@t: file type
	static const TabReaderPar& FileParams(eType t) { return TypeAttrs[int(t)].FileParam; }

	// Get Mutex type by file type
	static Mutex::eType MutexType(eType t) { return TypeAttrs[int(t)].MtxType; }

} fformat;

class TxtFile
	/*
	 * Basic class 'TxtFile' implements a fast buffered serial (stream) reading/writing text files
	 * (in Windows and Linux standart).
	 * Supports reading/writing zipped (gz) files.
	 * Optimised for huge files.
	 * Restriction: the default size of buffer, setting as NUMB_BLK * BasicBlockSize,
	 * should be bigger than the longest line in file. Otherwise file become invalid,
	 * and less than size_t.
	 * If size of reading files is less than default buffer's size,
	 * the buffer's size sets exactly to be sufficient to read a whole file.
	 * The reading/writing unit is a 'record', which in common case is a predifined set of lines.
	 * For FQ files 'record' is a set of 4 lines.
	 * For common text files 'record' is identical to 'line'.
	 * Empty lines are skipping by reading and therefore will not be writing by 'cloning' a file.
	 */
{
public:
	enum class eAction {
		READ,		// reads only existing file
		WRITE,		// creates file if it not exist and writes to it; file is cleared before
		READ_ANY	// creates file if it not exist and reads it
	};

protected:
	enum eFlag {			// signs of file
		// The first two right bits are reserved for storing the length of the LF marker: 1 or 2 
		// The first bit is set to 1 in the constructor.
		// If CR symbol is found at the end of line,
		//	the second bit is raised to 1, the first turn down to 0
		ISCR	  = 0x01,	// sign of Carriage Return ('\r')
		LFCHECKED = 0x02,	// the presence of a symbol CR is checked; for Reading mode
		ZIPPED	  = 0x04,	// file is zipped
		ABORTING  = 0x08,	// invalid file should be completed by throwing exception; for Reading mode
		ENDREAD	  = 0x10,	// last call of GetNextRecord() has returned NULL; for Reading mode
		PRNAME	  = 0x20,	// print file name in exception's message; for Reading mode
		MTHREAD	  = 0x40,	// file in multithread mode: needs to be locked while writing
		CLONE	  = 0x80,	// file is a clone
	};

	using bufflen = uint32_t;

private:
	constexpr static bufflen BlockSize = bufflen(8 * 1024 * 1024);	// 8 Mb
	const int CRLF = eFlag::ISCR + eFlag::LFCHECKED;	// combined 'flag'

	size_t	_fSize;			// the length of uncompressed file; for zipped file more than
							// 4'294'967'295 its unzipped length is unpredictable
	string	_fName;			// file's name
	mutable short _flag;	// bitwise storage for signs included in eFlag

protected:
	void* _stream;			// FILE* (for unzipped file) or gzFile (for zipped file)
	char* _buff;			// basic I/O (read/write) buffer
	bufflen	_buffLen;		// the length of the basic I/O buffer
	mutable bufflen _currRecPos;	// start position of the last readed/writed record in the block
	mutable size_t _recCnt;			// local counter of readed/writed records
	mutable Err::eCode	_errCode;
	//Stopwatch	_stopwatch;

	void SetFlag(eFlag f, bool val)	const { val ? _flag |= f : _flag &= ~f; }
	void RaiseFlag(eFlag f)	const { _flag |= f; }
	bool IsFlag(eFlag f)	const { return _flag & f; }
	bool IsZipped()			const { return _flag & ZIPPED; }
	bool IsClone()			const { return _flag & CLONE; }

	// *** 3 methods used by TxtReader only

	// Gets the number of characters corresponded to LF
	// In Windows LF matches '\r\n', in Linux LF matches '\n',
	//	return: in Windows always 1, in Linux: 2 for file created in Windows, 1 for file created in Linux
	BYTE LFSize() const { return static_cast<BYTE>(2) - (_flag & ISCR); }

	// Returns true if LF size is not defined
	//bool IsLFundef() const { return !(_flag & LFCHECKED); }

	// Establishes the presence of CR symbol at the end of line.
	//	@c: if c is CR then the second bit is raised to 1, the first turn down to 0,
	//	so the value return by LFSZ mask is 2, otherwise remains in state 1
	void SetLF(char c) { if (c == CR)	_flag |= CRLF; }

private:
	// Initializes instance variables, opens a file, sets a proper error code
	//	@param fName: valid full name of file
	//	@param mode: opening mode
	//	@param fStream: clonable file stream or NULL
	//	@returns: true is success, otherwise false
	bool SetBasic(const string& fName, eAction mode, void* fStream);

	// Allocates memory for the I/O buffer with checking
	//	@returns: true is success, otherwise false
	bool CreateIOBuff();

	//void operator = (const TxtFile&); //assignment prohibition

protected:
	// Constructs an instance: allocates I/O buffer, opens an assigned file
	//	@param fName: valid full name of assigned file
	//	@param mode: opening mode
	//	@param msgFName: true if file name should be printed in the exception's message
	//	@param abortInvalid: true if invalid instance shold be completed by throwing exception
	TxtFile(const string& fName, eAction mode, bool msgFName = true, bool abortInvalid = true);

#ifdef _MULTITHREAD
	// Constructs a clone of an existing instance
	// Clone is a copy of opened file with its own separate I/O buffer.
	// Used for multithreading file recording
	//	@param file: opened file which is cloned
	TxtFile(const TxtFile& file);
#endif

	// Destructs the instance: close assigned file, releases I/O buffer.
	~TxtFile();

	// Sets error code and throws exception if it is allowed.
	//	@param code: error code
	//	@param senderSpec: sender specificator
	void SetError(Err::eCode code, const string& senderSpec = strEmpty, const string& spec = strEmpty) const;

	// Returns true if instance is valid
	bool IsGood() const { return _errCode == Err::NONE; }
	// Returns true if instance is invalid
	bool IsBad()  const { return _errCode != Err::NONE; }

	// Returns current error code.
	//Err::eCode ErrCode() const	{ return _errCode; }

	// Returns conditional file name: name if it's printable, otherwise empty string.
	const string& CondFileName() const { return IsFlag(PRNAME) ? _fName : strEmpty; }

	// Returns size of uncompressed file,
	// or size of compressed file if its uncompressed length is less than UINT_MAX.
	// So it can be used for estimatation only!
	size_t Length() const { return _fSize; }

	// Returns number of readed/writed records.
	size_t RecordCount() const { return _recCnt; }

	// Throws exception
	//	@param msg: exception message
	void ThrowExcept(const string& msg) const { Err(msg, CondFileName()).Throw(); }

	// Throws exception
	//	@param code: error code
	void ThrowExcept(Err::eCode code) const { Err(code, CondFileName().c_str()).Throw(); }

public:
	// Returns file name
	const string& FileName() const { return _fName; }
};

// 'TxtReader' represents TxtFile for reading
//	Ñorrectly reads records even if the latter does not end with LF (CR,LF) character (s)
class TxtReader : public TxtFile
{
	reclen	_recLen = 0;			// the length of record with LF marker
	bufflen	_readedLen = 0;			// number of actually readed chars in block
	reclen* _linesLen = nullptr;	// array of lengths of lines in a record
	BYTE	_recLineCnt;			// number of lines in a record

	// Raises ENDREAD sign an return NULL
	//char* ReadingEnded()		  { RaiseFlag(ENDREAD); return NULL; }

	// Reads next block
	//	@param offset: shift of start reading position
	//	@returns: 0 if file is finished; -1 if unsuccess reading; otherwhise number of readed chars
	int ReadBlock(const bufflen offset);

	// Checks the actually read block to see if it exceeds the I/O buffer limit.
	// If block is complete (exceeds buffer limit), moves the misplaced remainder at the beginning if the buffer.
	//	@param currLinePos: current reading position
	//	@param blankLineCnt: number of empty lines read
	//	@returns true if block is complete
	bool CompleteBlock(bufflen currLinePos, bufflen blankLineCnt);

	// Establishes the presence of CR symbol at the end of line.
	void DefineLF();

	// Fills I/O buffer with 0, beginning from @offset position
	//void ClearBuff(size_t offset = 0) { memset(_buff + offset, 0, _buffLen - offset); }

	// Gets string containing file name and current record number.
	//const string RecordNumbToStr() const { return LineNumbToStr(0); }

	// Gets line number
	//	@param lineInd: index of line in a record
	size_t LineNumber(BYTE lineInd) const { return (_recCnt - 1) * _recLineCnt + lineInd + 1; }

protected:
	// Constructs an TxtReader instance: allocates buffers, opens an assigned file.
	//	@param fName: valid full name of assigned file
	//	@param mode: opening mode: READ or READ_ANY
	//	@param cntRecLines: number of lines in a record
	//	@param msgFName: true if file name should be printed in an exception's message
	//	@param abortInvalid: true if invalid instance should be completed by throwing exception
	TxtReader(const string& fName, eAction mode, BYTE cntRecLines, bool msgFName = true, bool abortInvalid = true);

	~TxtReader() { if (_linesLen)	delete[] _linesLen; }

	// Returns record without control
	char* RealRecord() const { return _buff + _currRecPos - _recLen; }

	// Gets current record.
	const char* Record() const { return IsFlag(ENDREAD) ? NULL : RealRecord(); }

	// Reads multiline record
	//	@returns: pointer to line or NULL if no more lines
	const char* GetNextRecord();

	// Reads one-line N-controlled record
	//	@param counterN: counter of 'N'
	//	@returns: pointer to line or NULL if no more lines
	const char* GetNextRecord(chrlen& counterN);

	// Reads one-line tab-controlled record
	//	@param tabPos: TAB's positions array that should be filled
	//	@param cntTabs: maximum number of TABS in TAB's positions array
	//	@returns: pointer to line or NULL if no more lines
	char* GetNextRecord(short* const tabPos, const BYTE tabCnt);

	// Returns line length in a multiline record
	//	@param lineInd: index of line in a record
	//	@param withoutLF: if true then without LF marker length
	reclen LineLengthByInd(BYTE lineInd, bool withoutLF = true) const {
		return _linesLen[lineInd] - withoutLF * LFSize();
	}

	// Throw exception if no record is readed.
	//void CheckGettingRecord() const {
	//	if (IsLFundef())
	//		ThrowExcept("attempt to get record's oinfo without reading record");
	//}

	// Returns the read pointer to the beginning of the last read record and decreases record counter
	//	Zeroes length of current reading record!
	//	@param sep: field separator character
	void RollBackRecord(char sep);

public:
	// Gets length of current reading record including LF marker
	//	Returns 0 after RollBackLastRecord() invoke.
	reclen RecordLength() const { return _recLen; }

	// Gets current readed line
	const char* Line()	const { return Record(); }

	// Returns string containing the file name and current line read number
	//	@param code: code of error occurs
	//	@param lineInd: index of line in a record; if 0, then first line
	const string LineNumbToStr(Err::eCode code = Err::EMPTY, BYTE lineInd = 0) const;

	// Throws exception occurred in the current reading line 
	//	@code: exception code
	//void ThrowLineExcept(Err::eCode code) const {
	//	Err(_errCode=code, LineNumbToStr().c_str()).Throw();
	//}

	// Throws exception with message included current reading line number
	//	@param msg: exception message
	void ThrowExceptWithLineNumb(const string& msg) const { Err(msg, LineNumbToStr().c_str()).Throw(); }

	// Gets length of current line without LF marker: only for single-line record! Used in Fa only??
	reclen LineLength()	const { return RecordLength() - LFSize(); }
};

#ifdef _TXT_WRITER
// 'TxtWriter' represents TxtFile for writing
class TxtWriter : public TxtFile
{
#ifdef _ISCHIP
	friend class ReadWriter;
#endif
public:
	static bool	Zipped;				// true if filed should be zippped

private:
	//typedef void(TxtWriter::* fAddChar)();
	//// 'Add delimiter' methods: [0] - empty method, [1] - adds delimiter and increases current position
	//static fAddChar fLineAddChar[2];

	char	_delim;					// field delimiter in line
	// === line write buffer
	char*	_lineBuff = NULL;		// line write buffer; for writing mode only
	reclen	_lineBuffLen = 0;		// length of line write buffer in writing mode, otherwise 0
	reclen	_lineBuffOffset = 0;	// current shift from the _buffLine; replaced by #define!!!
	BYTE	_floatBuffLen = 6;		// length of buffer where the resulting float as string is stored
	string	_formatFloat[2]{ "%4.2f", "%1.f" };	// format string for converting float to string
#ifdef _MULTITHREAD
	// === total counter of writed records
	size_t* _totalRecCnt;	// pointer to total counter of writed records; for clone only
	Mutex::eType _mtype;
#endif

	//void AddCharEmpty() {}
	////	Adds delimiter to the current position in the line write buffer and increases current position.
	//void AddDelim() { LineAddChar(_delim); }

	//	Adds delimiter
	//	@param addDelim: if true then adds delimiter to the current position in the line write buffer and increases current position
	void LineAddDelim(bool addDelim) { if (addDelim) LineAddChar(_delim); }

	// Allocates memory for write line buffer with checking
	//	@param len: size of the allocated buffer
	//	@returns: true if successful
	bool CreateLineBuff(reclen len);

	// Closes adding record to the I/O buffer: set current rec position and increases rec counter
	//	@param len: length of added record
	void EndRecordToIOBuff(bufflen len);

protected:
	// Returns current buffer write position
	reclen CurrBuffPos() const { return _lineBuffOffset; }

	// Adds character to the current position in the line write buffer
	//	@param ch: char to be set
	void LineAddChar(char ch) { _lineBuff[_lineBuffOffset++] = ch; }

	// Adds character to the current position in the line write buffer with adding delimiter
	//	@param ch: char to be set
	//	@param addDelim: if true then adds delimiter and increases current position
	void LineAddChar(char ch, bool addDelim);

	// Copies block of chars to the current position in the line write buffer,
	//	adds delimiter after string	and increases current position.
	//	@param src: pointer to the block of chars
	//	@param len: number of characters
	//	@param addDelim: if true then adds delimiter and increases current position
	//	@returns: new current position
	reclen LineAddChars(const char* src, reclen len, bool addDelim = true);

	// Constructs an instance: allocates I/O buffer, opens an assigned file.
	//	@param ftype: type of file
	//	@param fName: valid file name without extention
	//	@param printName: true if file name should be printed in the exception's message
	//	@param abortInvalid: true if invalid instance shold be completed by throwing exception
	TxtWriter(FT::eType ftype, const string& fName,
		char delim = TAB, bool printName = true, bool abortInvalid = true) :
		_delim(delim),
		_mtype(FT::MutexType(ftype)),
		TxtFile(fName + FT::Ext(ftype, Zipped), eAction::WRITE, printName, abortInvalid)
	{
#ifdef _MULTITHREAD
		_totalRecCnt = &_recCnt;	// for atomic increment
#endif
	}		// line buffer will be created in SetLineBuff()

#ifdef _MULTITHREAD
	// Constructs a clone of an existing instance (for multithreading file recording)
	// Clone is a copy of opened file with its own separate basic I/O and write line buffers
	//	@param file: opened file which is cloned
	TxtWriter(const TxtWriter& file);
#endif

	// Writes nonempty buffer, deletes line buffer and closes file
	~TxtWriter();

	// Sets current recording position of the line write buffer
	void LineSetOffset(reclen offset = 0) { _lineBuffOffset = offset; }

	// Increases current recording position
	//	@param len: size by which the current position is increased
	void LineIncrOffset(reclen len) { _lineBuffOffset += len; }

	template<typename... Args>
	void LineAddArgs(const char* format, Args ... args) {
		_lineBuffOffset += sprintf(_lineBuff + _lineBuffOffset, format, args...);
	}

	// Copies the string to the current position of the line write buffer
	//	@param str: string to be copied
	//	@param addDelim: if true then adds delimiter and increases current position
	//	@returns: new current position
	reclen LineAddStr(const string& str, bool addDelim = true) {
		return LineAddChars(str.c_str(), reclen(str.length()), addDelim);
	}

	// Adds two integral values separated by default delimiter to the current position of the line write buffer
	//	@param val1: first value to be set
	//	@param val2: second value to be set
	//	@param addDelim: if true then adds delimiter and increases current position
	void LineAddInts(ULONG val1, ULONG val2, bool addDelim = true);

	// Adds three unsigned integral values separated by default delimiter to the current position of the line buffer
	//	@param val1: first value to be set
	//	@param val2: second value to be set
	//	@param val3: third value to be set
	//	@param addDelim: if true then adds delimiter and increases current position
	void LineAddUInts(chrlen val1, chrlen val2, chrlen val3, bool addDelim = true);

	// Adds floating point value to the current position of the line write buffer,
	//	adds delimiter after value and increases current position.
	//	@val: value to be set
	//	@ndigit: number of digits to generate
	//	@addDelim: if true then adds delimiter and increases current position
	//void LineAddFloat(float val, BYTE ndigit, bool addDelim = true);

public:
	// Adds integral value to the current position of the line write buffer
	//	@param val: value to be set
	//	@param addDelim: if true then adds delimiter and increases current position
	void LineAddInt(LLONG val, bool addDelim = true);

	// Adds floating point value to the current position of the line write buffer
	//	@param val: value to be set
	//	@param addDelim: if true then adds delimiter and increases current position
	void LineAddFloat(float val, bool addDelim = true);

protected:
	// Adds single floating point value to the I/O buffer
	void LineAddSingleFloat(float val);

	// Adds string with a single value to IO buffer without checking buffer exceeding.
	//template<typename T>
	//void LineAddSingleValue(T val) { 
	//	LineAddStr(to_string(val), false);
	//	LineToIOBuff();
	//}

	// Adds line to the I/O buffer
	//	@param offset: start position for the next writing cycle
	void LineToIOBuff(reclen offset = 0);

	// Adds record to the I/O buffer.
	// Generates exception if writing is fall.
	//	@param src: record
	//	@param len: length of the record
	void RecordToIOBuff(const char* src, bufflen len);

	// Adds string to I/O buffer without checking buffer limit
	void StrToIOBuff(const string&& str);

	// Adds commented string to the I/O buffer
	//	@param str: recorded string to be presented as a comment
	void CommLineToIOBuff(const string& str) { StrToIOBuff("# " + str); }

	// Allocates line write buffer.
	//	@param len: length of buffer
	void SetLineBuff(reclen len);

	// Adds to line 2 int values and adds line to the IO buff.
	//void WriteLine(ULONG val1, ULONG val2);

	// Adds to line C-string with delimiter, and adds line to the IO buff.
	//void WriteLine(const char* str);

	// Adds to line string with delimiter and int value, and adds line to the IO buff.
	//void WriteLine(const string& str, int val);

	// Writes thread-safely current block to file
	void Write() const;

	// Adds content of another file (concatenates)
	//	return: true if successful
	//bool	AddFile(const string fName);

public:
	// Sets the number of digits in the fractional part when converting float to a string in buffer
	//	@param fractDigitsCnt: number of digits in the fractional part
	void SetFloatFractDigits(BYTE fractDigitsCnt);
};
#endif	// _TXT_WRITER

#ifndef _FQSTATN

// 'TabReader' represents TxtReader consisting of lines with tab-separated filds
class TabReader : public TxtReader
	/*
	 * Only number of fields defined in constructor is processed.
	 * Number of fields can be equal 1, in that case ordinary plain text is processed.
	 */
{
	FT::eType	_fType;
	reclen	_lineSpecLen = 0;		// length of the line specifier; for reading only
	short* _fieldPos = nullptr;		// array of start positions in _currLine for each field
	char*  _currLine = nullptr;		// current readed line; for reading only
	size_t _estLineCnt = vUNDEF;	// estimated number of lines

	// Checks if filed is valid, otherwise throws an exception
	//	@param fInd: field index
	//	@returns: true if field is valid
	bool IsFieldValid(BYTE fInd) const;

	// Initializes new instance
	//	@param type: file bioinfo type
	//	@param estLineCnt: if true then estimate count of lines
	void Init(FT::eType type, bool estLineCnt);

	// Frees allocated memory
	void Release() { if (_fieldPos) delete[] _fieldPos; _fieldPos = nullptr; }

	// Returns the read pointer to the beginning of the last read line and decreases line counter
	void RollBackLine() { RollBackRecord(TAB); }

protected:
	// Sets estimation of number of lines base obased on the line read, and rolls line back
	void SetEstLineCount();

	// Sets estimation of number of lines by type (for willge)
	//	@param wigType: some of willge type
	void SetEstLineCount(FT::eType wigType) {
		_estLineCnt = Length() / FT::FileParams(wigType).AvrLineLen;
	}

	// Reinitializes an instance of given type, correct estimated number of lines if it's predefined
	void ResetType(FT::eType type);

public:
	// Creates new instance for reading
	//	@param fName: name of file
	//	@param type: file bioinfo type
	//	@param msgFName: true if file name should be printed in exception's message
	//	@param abortInvalid: true if invalid instance should be completed by throwing exception
	TabReader(
		const string& fName,
		FT::eType type = FT::UNDEF,
		eAction	mode = eAction::READ,
		bool estLineCnt = true,
		bool msgFName = true,
		bool abortInvalid = true
	) : TxtReader(fName, mode, 1, msgFName, abortInvalid), _fType(type)
	{
		if (mode != eAction::WRITE && IsGood())
			Init(_fType, estLineCnt);
	}

#if defined _ISCHIP && defined  _MULTITHREAD
	// Creates a clone of TabReader class (for multithreading file recording)
	// Clone is a copy of opened file with its own buffer for writing only
	//  @param file: opened file which is copied
	TabReader(const TabReader& file) :
		_fType(file._fType),
		_lineSpecLen(file._lineSpecLen),
		_fieldPos(file._fieldPos),
		_currLine(NULL),
		TxtReader(file) {}
#endif

	~TabReader() { Release(); }

	// Returns a pointer to the substring defined by key
	//	@param str: null-terminated string to search the key
	//	@param key: string to search for
	//	@returns: a pointer to the null-terminated substring after key, or NULL if key does not appear in str
	static const char* KeyStr(const char* str, const string& key);

	// Checks definition or declaration line for key
	//	@param str: null-terminated string to search the key
	//	@param key: string to search for
	//	@returns: a pointer to the null-terminated substring followed after the key
	const char* CheckSpec(const char* str, const string& key);

	// Returns required int value with check
	//	@param str: null-terminated string to search the key
	//	@param key: string to search for
	//	@returns: key value, or throws an exception if key does not appear in str
	chrlen GetIntKey(const char* str, const string& key) { 
		return atoui(CheckSpec(str, key) + 1);
	}

	// Gets file bioinfo type
	FT::eType Type() const { return _fType; }

	// Gets count of readed lines
	size_t Count() const { return RecordCount(); }

	// Returns estimated number of lines
	size_t EstLineCount() const { return _estLineCnt; }

	// Skip commented lines and returns estimated number of uncommented lines
	//ULONG GetUncommLineCount();

	// Reads first line and set it as current.
	// Throws exception if invalid
	//	@cntLines: returned estimated count of lines.
	//	It works properly only if lines are sorted by ascending, f.e. in sorted bed-files.
	//	return: current line
	//const char*	GetFirstLine(ULONG& cntLines);

	// Reads next line and set it as current
	//	@param checkTabs: it true then check the number of incoming fields (tabs)
	//	@returns: read line
	const char* GetNextLine(bool checkTabs = true);

	// Gets current line
	char* GetLine() const { return _currLine; }

	// Returns null-terminated string by field's index from current record without check up (for FqReader)
	//	@param fInd: field index
	const char* StrField(BYTE fInd)	const { return _currLine + _fieldPos[fInd]; }

	// Returns null-terminated string by field's index from current record with check up (for FqReader)
	//	@param fInd: field index
	const char* StrFieldValid(BYTE fInd)	const {
		return IsFieldValid(fInd) ? StrField(fInd) : nullptr;
	}

	// Reads float by field's index from current line without check up
	//	@param fInd: field index
	float FloatField(BYTE fInd)	const { return float(atof(StrField(fInd))); }

	// Reads float by field's index from current line with check up
	//	@param fInd: field index
	float FloatFieldValid(BYTE fInd) const {
		return IsFieldValid(fInd) ? FloatField(fInd) : vUNDEF;
	}

	// Reads unsigned integer by field's index from current line without check up
	//	@param fInd: field index
	chrlen UIntField(BYTE fInd)	const { return atoui(StrField(fInd)); }

	// Reads integer by field's index from current line without check up
	//	@param fInd: field index
	//int IntField(BYTE fInd)	const { return atoi(StrField(fInd)); }

	// Reads integer by field's index from current line with check up
	//	@param fInd: field index
	//int IntFieldValid(BYTE fInd)	const {	return IsFieldValid(fInd) ? IntField(fInd) : vUNDEF; }

	// Reads unsigned long by field's index from current line without check up.
	//	@param fInd: field index
	//size_t ULongField(BYTE fInd)	const { return atoul(StrField(fInd)); }

	// Reads long by field's index from current line without check up.
	//	@param fInd: field index
	//long LongField(BYTE fInd)	const { return atoul(StrField(fInd)); }

	// Reads long by field's index from current line with check up.
	//	@param fInd: field index
	//long LongFieldValid(BYTE fInd)	const {	return IsFieldValid(fInd) ? LongField(fInd) : vUNDEF; }

	// Initializes region by two consecutive integer fields
	//	@param fInd: first field index
	//	@param rgn: region that is initialized
	void InitRegion(BYTE fInd, Region& rgn) const;
};


// 'ChromDefRegions' represents a container of defined regions within chromosome
// Defined in TxtFile.h since it's used in Fa class
class ChromDefRegions : public Regions
{
	// 'Combiner' consistently combines regions separated by a space less than the minimum allowed
	class Combiner
	{
		chrlen	_gapLen;	// minimum allowable gap between regions
		Region	_rgn;		// current region

	public:
		// Constructor accepting the minimum allowable gap
		//	@param minGapLen: minimal length which defines gap as a real gap
		Combiner(chrlen minGapLen) : _gapLen(minGapLen) {}

		// Returns false if the next region is separated from the current one
		// by a value less than the established minimum gap.
		// Otherwise returns true and combined region as parameter.
		bool ExceptRegion(Region& rgn);

		// Returns last excepted region
		const Region& LastRegion() const { return _rgn; }
	};

	const size_t DefCapacuty = 12;	// default container capacity

	string	_fName;
	chrlen	_gapLen;	// total length of gaps
	mutable bool _new;

public:
	static const string Ext;	// regions file extension

	// Creates new instance and initializes it from file if one exist
	//	@param fName: chrom file name without extension
	//	@param minGapLen: length, gaps less than which are ignored when reading; if 0 then read all regions 
	ChromDefRegions(const string& fName, chrlen minGapLen = 2);

	// Returns true if instance is empty
	bool Empty() const { return _new; }

	// Returns total length of gaps
	chrlen GapLen() const { return _gapLen; }

	// Increments total length of gaps
	void IncrGapLen(chrlen gapLen) { _gapLen += gapLen; }

	// Adds defined region beginning from current gap position
	//	@param rgn: defined region
	//	@param minGapLen: minimal length which defines gap as a real gap
	void AddRegion(const Region& rgn, chrlen minGapLen);

	// Saves new instance
	void Write() const;

#if defined _READDENS || defined _BIOCC
	// Combines regions with a gap less than a specified value
	//	@param minGapLen: minimal length which defines gap as a real gap
	void Combine(chrlen minGapLen);
#endif
};

// 'FaReader' supports reading chromosomes in FA format
class FaReader : public TxtReader
{
	// 'DefRgnMaker' produced chrom defined regions while FaReader reading 
	class DefRgnMaker
	{
		chrlen	_minGapLen,			// minimal length which defines gap as a real gap
				_currPos;			// current readed nt position
		Region	_defRgn;			// current defined region
		ChromDefRegions& _defRgns;	// external chrom's defined regions

	public:
		DefRgnMaker(ChromDefRegions& rgns, chrlen minGapLen)
			: _minGapLen(minGapLen), _currPos(0), _defRgns(rgns) { _defRgns.Clear(); }

		// Adds last readed line length
		void AddLineLen(chrlen lineLen) { _currPos += lineLen; }

		// Remove last readed line length
		void RemoveLineLen(chrlen lineLen) { _currPos -= lineLen; }

		// Adds a gap to assign defined regions
		//	@param start: gap's start position.
		//	@param len: gap's length
		void AddGap(chrlen start, chrlen len);

		// Closes adding gaps, saves hrom's defined regions
		//	@param cLen: chrom length
		void CloseAddGaps(chrlen cLen);
	};

	const char FaComment = '>';
	chrlen		_cLen;						// length of chromosome
	const char* (FaReader::* _pGetLine)();	// pointer to the 'GetNextLine' method
	unique_ptr<DefRgnMaker> _rgnMaker;		// chrom defined regions store

	// Searches 'N' subsequence in the current line beginning from the start position
	// and adds complete 'N' subsequence to _rgnMaker
	//	@param startPos: starting line reading position
	//	@param NCnt: number of 'N' in the line beginning from the start position
	//	@returns: true if line trimming is complete, false for single 'N'
	void CountN(chrlen startPos, reclen NCnt);

	// Reads line and set it as current with def filling regions
	// Since method scans the line, it takes no overhead expenses to do this.
	const char* GetLineWitnNControl();

public:
	// Opens existing file, reads the first line and sets chrom length
	//	@param fName: file name
	//	@param rgns: def regions to fill, otherwise nullptr to reading without 'N' control
	FaReader(const string& fName, ChromDefRegions* rgns = nullptr);

	// Gets chromosome's length
	chrlen ChromLength() const { return _cLen; }

	// Reads line and set it as current (with filling def regions if necessary)
	const char* NextGetLine() { return (this->*_pGetLine)(); }

	// Closes reading
	void CLoseReading() { if (_rgnMaker) _rgnMaker->CloseAddGaps(_cLen); }
};

#endif	// no _FQSTATN
