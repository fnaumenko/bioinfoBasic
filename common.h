/**********************************************************
common.h 
Provides common functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 07/28/2024
***********************************************************/
#pragma once

#ifndef _COMMON_H
#define _COMMON_H

#include <stdint.h>		// uint32_t
#include <string>
#include <iostream>	
#include <sstream>		// for NSTR( x )
#include <iomanip>		// setprecision(), setw()
#include <vector>
#include <sys/stat.h>	// struct stat
#include <limits>       // std::numeric_limits
#include <functional>
#include <chrono>
#ifdef _MULTITHREAD
#include <mutex>
#endif

#ifdef __unix__
	#include <unistd.h>
	#include <stdlib.h>
	#include <memory>
	#include <limits.h>
	//#include <math.h>
	#include <dirent.h>
	#include <stdio.h>
#ifdef _MULTITHREAD
	//#include <pthread.h>
	#define InterlockedExchangeAdd	__sync_fetch_and_add
	#define InterlockedIncrement(p)	__sync_add_and_fetch(p, 1)
#endif
	#include <string.h>		// strerror_r()
	#include <stdexcept>	// throw std exceptions

	typedef unsigned char	BYTE;
	typedef	unsigned short	USHORT;
	typedef	unsigned int	UINT;
	typedef unsigned long	ULONG;
	typedef unsigned long long ULLONG;
	typedef long long LLONG;
	typedef struct stat struct_stat64;
	//typedef void*	thrRetValType;

	#define _stricmp strcasecmp	// case-sensitive comparison
	#define _fseeki64 fseeko64
	#define _ftelli64 ftello64
	#define _fileno	fileno
	#define _gcvt	gcvt
	#define _stat64 stat
	#define fopen	fopen64
	#define TRUE	1
	#define FALSE	0
	#define isNaN(x)	(x!=x)	// will not work by using --ffast-math compiler flag
	#define	LOCALE_ENG	"en_US.utf8"
	//#define SLASH '/'	// standard Linux path separator
	//#define F_READ_MODE O_RDONLY
#elif defined _WIN32
	#define OS_Windows
	#include <windows.h>
//#ifdef _MULTITHREAD
	//#define pthread_mutex_t CRITICAL_SECTION
//#endif	
	typedef unsigned __int64 ULLONG;
	typedef __int64 LLONG;
	typedef struct __stat64 struct_stat64;

	#define isNaN isnan
	#define	LOCALE_ENG	"English"
#endif		// OS Windows

#ifdef _ZLIB
	#ifdef OS_Windows
		#define ZLIB_WINAPI
	#endif
	#include "zlib.h"
	#if ZLIB_VERNUM >= 0x1240
		#define ZLIB_NEW
	#else	
		#define ZLIB_OLD
	#endif
#endif	// _ZLIB

// specific types
typedef BYTE		thrid;		// type number of thread
typedef BYTE		chrid;		// type number of chromosome
typedef uint16_t	readlen;	// type length of Read
typedef uint32_t	chrlen;		// type length of chromosome
typedef chrlen		fraglen;	// type length of fragment
typedef chrlen		genlen;		// type length of genome

#define	CHRLEN_UNDEF	-1			// undefined length of chromosome
#define	CHRLEN_MAX		LONG_MAX //-1	// max length of chromosome
//#define CHRLEN_CAPAC	10	// capacity of max chrom length;
							// may be count by DigitsCountUInt() every time,
							// but is the same if chrlen defined as int or long,
							// so is defined as static value

#define cNULL	'\0'
#define cN		'N'
#define HPH		'-'
#define USCORE	'_'
#define SPACE	' '
#define sSPACE	" "
#define QUOT	'\''
#define DOT		'.'
#define COLON	':'
#define COMMA	','
#define PERS	'%'
#define AT		'@'
#define PLUS	'+'
#define HASH	'#'
#define TAB		'\t'
#define LF		'\n'	// Line Feed, 10
#define CR		'\r'	// Carriage Return, 13
//#define sBACK	"\b"

using namespace std;

static const chrlen CHRLEN_CAPAC = numeric_limits<chrlen>::digits10 + 1;

static const char* SepCl = ": ";		// colon separator
static const char* SepSCl = "; ";		// semicolon separator
static const char* SepCm = ", ";		// comma separator
static const char* SepDCl = ":: ";		// double colon separator

static const std::string ZipFileExt = ".gz";
static const string strEmpty = "";

/*** COMMON MACROS & FUNCTION ***/

#define vUNDEF	-1	// undefined value
#define NO_VAL	-1	// option value is prohibited
#define NO_DEF	-1	// do not print option default value


// Gets number of members in static array
#define ArrCnt(arr)	sizeof(arr)/sizeof(arr[0])

static const char* Booleans[] = { "OFF","ON" };	// boolean values
// Returns C-string 'ON' or 'OFF'
inline const char* BoolToStr(bool val) { return Booleans[val]; }

//===== operations with integers

// Returns number of ones in the bit representation of an integer
int OnesCount(int n);

// Returns right position of right one in the bit representation of an integer
int RightOnePos(int n);

// Gets number of digist in a unsigned integral value
//	@param val: unsigned integral value
//	@returns: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCountUInt (uint32_t val);

// Gets number of digist in a integral value
//	@param val: integral value
//	@returns: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCountInt(int32_t val);

// Gets number of digist in a long integral value
//	@param val: long integral value
//	@returns: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCountULong(uint64_t val);

// Gets number of digist in a integral value
//	@param val: integral value
//	@param isLocale: if true then adds number of '1000' separators
//	@returns: number of digist without minus symbol or 0 if value is 0
BYTE DigitsCountLocale(uint64_t val, bool isLocale = false);


// Converts string to chrlen
//	@param str: reference to C-string beginning with the representation of an unsigned integral number.
//	Must start with a digit.
chrlen atoui_by_ref(const char*& str);

// Converts string to chrlen
//	@param str: C-string beginning with the representation of an unsigned integral number.
//	Must start with a digit.
inline chrlen atoui(const char* str) { return atoui_by_ref(str); }

// Converts string to size_t
//	@param str: C-string beginning with the representation of an unsigned long integral number.
//	Must start with a digit.
size_t atoul(const char* str);

// Digital to STRing
// Returns value's string representation. http://rootdirectory.de/wiki/NSTR()
// Instead of std::to_string(x) [C++11] because of compatibility
//#define NSTR(x) std::to_string(x)

// ===== miscellaneous

// Returns percent of part value relatively total value
//	@param part: part of the total value
//	@param total: total value
inline float Percent(size_t part, size_t total) { return total ? 100.f * part / total : 0.f; }

// Returns string represents the percent of part relatively total
//	@param percent: value of percent
//	@param precision: count of fractional digits; 
//	if count of value's mapped digits is more then that, printed "<X%", or exactly by default
//	@param fieldWith: displayed width of value and '%' or '<' sign (excluding parentheses), or exactly if 0;
//	@param parentheses: if true then parenthesize the value (not considering fieldWith)
string PercentToStr(float percent, BYTE precision = 0, BYTE fieldWith = 0, bool parentheses = false);

// Returns string represents the percent of part relatively total
//	@param part: value represents desired %
//	@param total: value represents 100%
//	@param precision: count of mapped digits; 
//	if count of value's mapped digits is more then that (too little percent), printed "<n%" or exactly by default
//	@param fieldWith: the width of the display field insine parentheses or exactly by default; should include a '%' mark
//	@param parentheses: if true parenthesize the value
inline string sPercent(size_t part, size_t total, BYTE precision=0, BYTE fieldWith=0, bool parentheses=false) {
	return PercentToStr(Percent(part, total), precision, fieldWith, parentheses);
}

// Gets linear density (density per 1000 bs) of some elements
//	@param cnt: number of elements
//	@param len: length on which the density is determined
//	Used in isChIP (Imitator.cpp) and readDens (readDens.h).
inline float LinearDens(size_t cnt, chrlen len) { return len ? 1000.f * cnt / len : 0; }

// Prints solid horizontal line
//	@param lw: width of line
void PrintSolidLine(USHORT lw);

#ifdef _BIOCC

// Align position to the up or down resolution level
// f.e. by resoluation==5 pos 102 -> 100+relative, pos 104 -> 105++relative
//	@param pos: chromosome's position
//	@param res: resolution
//	@param relative: 0 (0-relative position, BED) or 1 (1-relative position: the first base is 1, WIG)
//	@returns: aligned position
chrlen AlignPos(chrlen pos, BYTE res, BYTE relative);

#endif

//#ifdef OS_Windows
//string	Wchar_tToString(const wchar_t* wchar);
//wchar_t* StringToWchar_t(const string &str, wchar_t* wchar);
//#endif

/*** end of COMMON FUNCTION ***/

#ifdef _DUP_OUTPUT
#include <fstream>
// 'dostream' duplicates outstream to stdout & file
class dostream : public std::ostream
{
	std::ofstream file;
public:
	//dostream(std::ostream& s) : std::ostream(cout.rdbuf()) {}

	dostream() : std::ostream(cout.rdbuf()) {}

	~dostream() { if (file.is_open())	file.close(); }		// in case of exception or holding execution by user

	std::ofstream& File() { return file; }

	// Open output file with given name
	//	@returns: true if the attempt was unsuccessful, false otherwise
	bool OpenFile(const string& fname);

	// Associates 'loc' to the dostream as the new locale object to be used with locale-sensitive operations
	void Imbue(const std::locale& loc);

	template <typename T> dostream& operator<< (T val) {
		cout << val;
		file << val;
		return *this;
	}
};

extern dostream dout;		// stream's duplicator
#else
	#define dout	cout
	#define dostream ostream
#endif	// _DUP_OUTPUT

// 'Product' keeps Product oinfo
static struct Product
{
	static const string	Title;
	static const string	Version;
	static const string	Descr;

	// Gets length of product name (title plus version).
	//static BYTE MarkLength() { return BYTE(Title.length()); // + strlen(Version) + 1); }

	// Gets title plus version.
	//static const string& Name() { return Title + string(1, HPH) + string(Version); }
} product;


typedef pair<float,float> pairVal;	// pair of float values

	// 'PairVals' holds pair of values and their min and max limits
class PairVals
{
	pairVal vals[3];

public:
	enum Type { SET = 0, MIN = 1, MAX = 2 };

	PairVals(float val1, float val2, float min1, float min2, float max1, float max2);

	// Gets pair of values
	const pairVal& Values(Type t = SET) const { return vals[t]; }
};

#define ErrWARNING	Err::NONE

// 'Error' implements error's (warning's) messages and treatment.
class Err
{
public:
	enum eCode {
		NONE,		// warning message
		MISSED,		// something missed
		F_NONE,		// no file
		D_NONE,		// no directory
		FD_NONE,	// nor file nor directory 
		F_MEM,		// TxtFile(): memory exceeded
		F_OPEN,		// TxtFile(): file open error
		F_CLOSE,	// TxtFile(): file close error
		F_READ,		// TxtFile(): file read error
		F_EMPTY,
		F_BIGLINE,	// TxtFile(): too big line
		FZ_MEM,		// TxtFile(): not enough gzip buffer
		FZ_OPEN,	// TxtFile(): gzip open error
		FZ_BUILD,	// TxtFile(): doen't support gzip
		F_WRITE,	// file write error
		//F_FORMAT,
#ifndef _FQSTATN
		TF_FIELD,	// TabReader: number of fields is less than expected
		TF_EMPTY,	// TabReader: no records
#endif
		EMPTY
	};

private:
	static const char* _msgs[];
	enum eCode	_code;			// error code
	char * _outText;			// output message
	//unique_ptr<char[]> _outText;			// output message

	// Initializes _outText by C string contained message kind of
	// "<sender>: <text> <specifyText>".
	void set_message(const char* sender, const char* text, const char* specifyText=NULL);

	//void set_message(const string& sender, const char* text, const char* specifyText=NULL) {
	//	set_message(sender==strEmpty ? NULL : sender.c_str(), text, specifyText);
	//}

public:
	// Returns string containing file name and issue number.
	//	@issName: name of issue
	//	@issNumb: number of issue
	//	@fName: file name
	//static const string IssueNumbToStr(const string& issName, ULONG issNumb, const string& fName);

	// Gets message "no <filename>.ext[.gz] file(s) in this directory"
	//	@param fName: file name or template
	//	@param plural: if true then in plural form
	//	@param fExt: file extention or none
	static const string MsgNoFile(const string & fName, bool plural, const string fExt = strEmpty);

	// Code-attached constructor.
	//	@param code: exception/warning message as code
	//	@param sender: name of object who has generated exception/warning, or NULL if no sender
	//	@param specifyText: aditional text to specify exception/warning message
	Err(eCode code, const char* sender, const char* specifyText=NULL): _code(code) {
		set_message(sender, _msgs[code], specifyText);
	}

	// Code-attached constructor.
	//	@param code: exception/warning message as code
	//	@param sender: name of object who has generated exception/warning, or NULL if no sender
	//	@param specifyText: aditional text to specify exception/warning message
	Err(eCode code, const char* sender, const string& specifyText) : _code(code) {
		set_message(sender, _msgs[code], specifyText.c_str());
	}

	// C-string-attached constructor.
	//	@param text: exception/warning message
	//	@param sender: name of object who has generated exception/warning
	Err(const char* text, const char* sender=NULL) : _code(NONE) {
		set_message(sender, text);
	}

	// String-attached constructor.
	//	@param text: exception/warning message
	//	@param sender: name of object who has generated exception/warning
	Err(const char* text, const string& sender) : _code(NONE) {
		set_message(sender.c_str(), text);
	}

	// String-attached constructor.
	//	@param text: exception/warning message
	//	@param sender: name of object who has generated exception/warning
	Err(const string& text, const char* sender=NULL) : _code(NONE) {
		set_message(sender, text.c_str());
	}

	// String-attached constructor.
	//	@param text: exception/warning message
	//	@param sender: name of object who has generated exception/warning
	Err(const string& text, const string& sender) : _code(NONE) {
		set_message(sender.c_str(), text.c_str());
	}

	// copy constructor
	Err(const Err& src);

	~Err() { delete[] _outText; _outText = NULL; }

	const char* what() const /*throw()*/ { return _outText; }

	eCode Code() const		{ return _code; }

	// Returns point to the aditional text to specify exception/warning message or NULL
	//const char* SpecifyText() const	{ return _specifyText; }

	// Throws exception or outputs Err message.
	//	@param throwExc: if true then throws exception, otherwise outputs Err message
	//	@param eol: if true then carriage return after Err message
	void Throw(bool throwExc = true, bool eol = true);
	
	// Outputs warning
	//	@param prefix: output ": " before "WARNING"
	//	@param eol: if true then carriage return after Err message
	void Warning(bool eol = true, bool prefix = false);
};

// 'FileSystem' implements common file system routines
static class FS
{
	// Returns true if file system's object exists
	//	@param name: object's name
	//	@param st_mode: object's system mode
	static bool IsExist(const char* name, int st_mode);

	// Checks if file system's object doesn't exist
	//	@param name: object's name
	//	@param st_mode: object's system mode
	//	@param throwExcept: if true then throws exception, otherwise prints warning without LF
	//	@param ecode: error's code
	//	@returns: true if file or directory doesn't exist (make sence while 'throwExcept' set to false)
	//	@throws 'ecode' if system's object doesn't exist
	static bool CheckExist	(const char* name, int st_mode, bool throwExcept, Err::eCode ecode);

	// Searches through a file name for the any extention ('/' or '\' insensible).
	//	@param name: name of file
	//	@returns: index of the DOT matched extention; otherwise string::npos
	static size_t GetLastExtPos	(const string &name);

	// Returns true if file name has specified extention ignoring zip extention (case-insensitive)
	//	@param name: name of file
	//	@param ext: file extention
	//	@param knownZip: if true then file name has ZIP extention
	//	@param composite: if true then file has a compound extention
	//	@returns: true if file name has specified extention
	static bool HasCaseInsExt(const string &name, const string &ext, bool knownZip, bool composite = true);

	// Returns file name without extention ('/' or '\' insensible)
	static string const FileNameWithoutLastExt (const string& name) {
		return name.substr(0, GetLastExtPos(name));
	}

public:
	// === file size

	// Returns size of file or -1 if file doesn't exist
	static LLONG Size 	(const char* name);

	// Returns real size of zipped file or -1 if file cannot open; limited by UINT
	static LLONG UncomressSize	(const char* name);

	// === check dir/file existing

	// Returns true if file exists
	static bool IsFileExist	 (const char* name) { return IsExist(name, S_IFREG); }
	
	// Returns true if directory exists
	static bool IsDirExist	 (const char* name) { return IsExist(name, S_IFDIR); }
	
	// Returns true if file or directory exists
	static bool IsFileDirExist(const char* name) { return IsExist(name, S_IFDIR|S_IFREG); }
	
	// Checks if file doesn't exist
	//	@param name: name of file
	//	@param throwExcept: if true then throws exception, otherwise prints warning without LF
	//	@returns: true if file doesn't exist (make sence while 'throwExcept' set to false)
	//	@throws Err::F_NONE if file doesn't exist
	static bool CheckFileExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFREG, throwExcept, Err::F_NONE);
	}

	// Checks if directory doesn't exist
	//	@param name: name of directory
	//	@param throwExcept: if true then throws exception, otherwise prints warning without LF
	//	@returns: true if directory doesn't exist (make sence while 'throwExcept' set to false)
	//	@throws Err::D_NONE if directory doesn't exist
	static bool CheckDirExist	(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR, throwExcept, Err::D_NONE);
	}

	// Ñhecks file or directory for existence
	//	@param name: name of file or directory
	//	@param throwExcept: if true then throws exception, otherwise prints warning without LF
	//	@returns: true if file or directory doesn't exist (make sence while 'throwExcept' set to false)
	//	@throws Err::FD_NONE if file or directory doesn't exist
	static bool CheckFileDirExist(const char* name, bool throwExcept = true) {
		return CheckExist(name, S_IFDIR|S_IFREG, throwExcept, Err::FD_NONE);
	}

	// Ñhecks file or directory for existence
	//	@param name: name of directory
	//	@param ext: file extention; if set, check for file first
	//	@param throwExcept: if true then throws exception, otherwise prints warning without LF
	//	@returns: true if file or directory doesn't exist (make sence while 'throwExcept' set to false)
	//	@throws Err::F_NONE or Err::D_NONE if file or directory doesn't exist
	static bool CheckFileDirExist(const char* name, const string & ext, bool throwExcept);

	// === check dir/file name

	// Returns the name of an existing file or directory
	//	@param name: name of file or directory
	//	@returns: input parameter 'name'
	//	@throws Err::FD_NONE if file or directory does not exist
	static const char* CheckedFileDirName(const char* name);

	// Returns the name of an existing file
	//	@param name: name of file
	//	@returns: input parameter 'name'
	//	@throws Err::F_NONE if file does not exist
	static const char* CheckedFileName(const char* name);

	// Returns a pointer to the path checked if it exist, otherwise throws exception
	//	@param opt: Options value
	//	@returns:  pointer to the checked path
	//static const char* CheckedDirName	(int opt);

	// Returns true if directory is open to writing
	static bool IsDirWritable(const char* name) { 
#ifdef OS_Windows
		return true;
#else
		return !access(name, W_OK); 
#endif
	}

	// === check file extension

	// Returns true if file has any extension.
	//	@param fname: file name
	static bool HasExt	(const string &fname) { return GetLastExtPos(fname) != string::npos; }

	// Returns true if file has a specified  extension. Case insensitive search
	//	@param ext: extension includes dot symbol
	//	@param composite: true if the extension is real (i.e. strictly last)
	static bool HasExt	(const string& fname, const string& ext, bool composite = true) {
		return HasCaseInsExt(fname, ext, HasGzipExt(fname), composite);
	}

	// Returns true if file has '.gz' extension. Case insensitive search
	static bool HasGzipExt(const string& fname) { return HasCaseInsExt(fname, ZipFileExt, false); }

	// Returns string containing real file extension (without gzip).
	//	@param fname: pointer to the file name
	//	@returns:  string containing real file extension or empty string if no real extention
	static string const GetExt(const char* fname);

	// === dir/file name transform

	// Returns file name without extention ('/' or '\' insensible)
	static string const FileNameWithoutExt (const string& fname) {
		return FileNameWithoutLastExt( HasGzipExt(fname) ? FileNameWithoutLastExt(fname) : fname );
	}

	// Returns true if file name is short (without path)
	static bool IsShortFileName(const string& fname);

	// Returns short file name by long one
	//	@param fname: long file name
	static string const ShortFileName (const string& fname);

	// Returns directory name by long file name
	//	@param fname: long file name
	//	@addSlash: true if slash sould be added at the end
	static string const DirName (const string& fname, bool addSlash = false);

	// Returns the name of last subdirectory by long file name
	//	@param name: long dir name
	static string const LastSubDirName (const string& name);

	// Returns the name of last subdirectory
	//	@param name: long dir name
	static string const LastDirName (const string& fname);

	// Returns the name ended by slash. 
	// If the name already ends with a slash, returns it.
	//	@param name: short or long name
	static string const MakePath(const string& name);

	// Returns composed file name
	//	@param oName: used as a base output short or long file name, or iutput directory name, or NULL. 
	//	If NULL then the input file name is used as a base 
	//	@param iName: input short or long file name
	//	@param suffix: output file suffix (can be extention, started with DOT), empty by default. 
	//	If output file extention matches the extension of the input file, the suffix '_out' is added to the resulting file name
	//	@throws Err::D_NONE if directory in name doesn't exist
	static string const ComposeFileName(const char* oName, const char* iName, const string& suffix = strEmpty);

	// === files in dir

#ifndef _FQSTATN
	// Fills external vector of strings by file's names found in given directory
	// Implementation depends of OS.
	//	@param files: external vector of strings that should be filled by file's names
	//	@param dirName: name of directory
	//	@param ext: file's extention as a choosing filter
	//	@param all: true if all files with given extention should be placed into external vector, otherwise only one (any)
	//	@returns:  true if files with given extention are found
	static bool GetFiles (vector<string>& files, const string& dirName, const string& ext, bool all = true);
#endif

//	static void	Delete		(const char* fname) {
//#ifdef OS_Windows
//		DeleteFile(fname);
//#else
//		unlink(fname);
//#endif
//	}
} fs;

// Basic class for wall time measurement
class TimerBasic
{
protected:
	// Prints elapsed wall time interval
	//	@param elapsed: elapsed time in milliseconds
	//	@param title: string printed before time output
	//	@param parentheses: if true then wrap the time output in parentheses
	//	@param prSec: if true and parentheses are not set then print "sec" after seconds with fractional part
	static void Print(long elapsed, const char *title, bool parentheses, bool prSec);

	mutable chrono::steady_clock::time_point	_begin;
	mutable bool	_enabled;	// True if local timing is enabled

	// Creates a new TimerBasic
	//	@param enabled: if true then set according total timing enabling
	TimerBasic(bool enabled = true) { _enabled = enabled ? Enabled : false;	}
	
	// Stops timer and return elapsed wall time in seconds
	long GetElapsed() const;
	
public:
	// True if total timing is enabled
	static bool		Enabled;

	// True if instance timing is enabled
	bool IsEnabled() const { return _enabled; }

	// Starts timer
	void Start() { if(_enabled)	_begin = chrono::steady_clock::now(); }
};

// 'Timer' measures the single wall time interval
class Timer : public TimerBasic
{
private:
	static clock_t	_StartCPUClock;

public:
	// Starts enabled CPU timer, if it is enabled
	static void StartCPU()	{ if( Enabled ) _StartCPUClock = clock(); }
	
	// Stops enabled CPU timer and print elapsed time
	//	@param isLF: if true then ended output by LF
	static void StopCPU(bool isLF=true) {
		if(Enabled)	Print( 1000 * (clock()-_StartCPUClock) / CLOCKS_PER_SEC, "CPU: ", false, isLF);
	}

	// Creates a new Timer and starts it if timing is enabled
	//	@param enabled: if true then set according total timing enabling
	Timer(bool enabled = true) : TimerBasic(enabled) { Start(); }
	
	// Stops enabled timer and prints elapsed time with title
	//	@param title: string printed before time output
	//	@param parentheses: if true then wrap the time output in parentheses
	//	@param prSec: if true and parentheses are not set then print "sec" after seconds with fractional part
	void Stop(const char *title, bool parentheses = false, bool prSec = true) {
		if(_enabled)	Print(GetElapsed(), title, parentheses, prSec);
	}

	// Stops enabled timer and prints elapsed time
	//	@param offset: space before time output
	//	@param parentheses: if true then wrap the time output in parentheses
	//	@param prSec: if true and parentheses are not set then print "sec" after seconds with fractional part
	void Stop(int offset = 0, bool parentheses = false, bool prSec = true);
};

#ifdef _TEST

// 'Stopwatch' measures the sum of wall time intervals
class Stopwatch : public TimerBasic
{
	private:
		mutable long	_sumTime;
		mutable bool	_isStarted;		// true if Start() was called even ones

	public:
		Stopwatch() : _sumTime(0), _isStarted(false), TimerBasic() {}

		// Starts Stopwatch
		void Start() const		{ ((TimerBasic*)this)->Start(); _isStarted = true; }

		// Stops Stopwatch
		//	@param title: if not empty, and if instance was launched, output sum wall time with title
		//	'const' to apply to constant objects
		void Stop(const string& title = strEmpty) const;
};

// 'Stopwatch' measures the sum of CPU time (clocks) intervals
class StopwatchCPU
{
	private:
		clock_t	_clock;
		clock_t	_sumclock;

	public:
		StopwatchCPU() : _sumclock(_clock=0) {}

		void Start(bool reset=false)	{ _clock = clock(); if(reset) _sumclock = 0; }

		// Stops StopwatchCPU
		//	@param title: string printed before time output
		//	@param print: if true time should be printed
		//	@param isLF: if true then ended output by LF
		void Stop(const char* title, bool print = false, bool isLF = false);
};

#endif	// _TEST

static class Mutex
{
public:
	enum class eType { OUTPUT, INCR_SUM, WR_BED, WR_SAM, WR_FQ, NONE };
#ifdef _MULTITHREAD
private:
	static bool	_active;	// true if multithreading is set
	//static pthread_mutex_t	_mutexes[];
	static mutex	_mutexes[];
public:
	// Returns true if Mutex type is not NONE
	static bool IsReal(eType mtype) { return mtype != eType::NONE;	}
	static void Init(bool active)	{ _active = active; }
	static bool isOn()				{ return _active; }
	static void Lock(eType type)	{ if (_active) _mutexes[int(type)].lock();	}
	static void Unlock(eType type)	{ if (_active) _mutexes[int(type)].unlock(); }
#endif
} myMutex;

// 'Chrom' establishes correspondence between chromosome's ID and it's name.
static class Chrom
/**********************************************************************************

===== TERMS ======
MARK: human view of chrom's number:	1, ..., Y
ABBREVIATION CHROM's NAME:			chr1, ..., chrY
SHORT CHROM's NAME:					chrom 1, ..., chrom Y
TITLE CHROM's NAME:					chromosome X, ..., chromosome Y
PREFIX: substring before mark:		chr | chrom | chromosome.

===== CHROM's ID NUMBERING DISCIPLINES =====
The autosomes' ID corresponds to the chrom number reduced by 1 (0-based numbering).
For the heterosomes' 2 different disciplines are applied:
RELATIVE: the heterosomes' ID starts with the last autosomes' ID increased by 1, in the next order: X, Y, M
ABSOLUTE: the heterosomes' ID is a code of mark character.

Relative discipline is applied when the number of the last autosome is known, i.e. when using chrom sizes data.
Otherwise it is impossible to determine the sequential ID of the first heterosome;
therefore, the code of heterosome's letter is used as its ID.

Absolute discipline (default) applies in case of sequential processing input BED.
Specifying chrom size data or processing input BAM automatically changes discipline to relative.

===== CHROM DATA DEPENDENCIES =====
User selectable chromosome (customization)		needs:	absolute | relative chrom ID discipline
Control of the item belonging to the chromosome	needs:	absolute | relative chrom ID discipline
Chrom size control								needs:	absolute | relative chrom ID discipline + ChromSizes
Using BAM										needs:	relative chrom ID discipline

Chrom class distinguishes these disciplines by the value of the 'firstHeteroID' variable:
a zero (default) value indicates an absolute discipline, a non-zero value indicates a relative one.

***********************************************************************************/
{
public:
	static const char*	Abbr;				// Chromosome abbreviation
	static const BYTE	MaxMarkLength = 2;	// Maximal length of chrom's mark
	static const BYTE	MaxAbbrNameLength;	// Maximal length of abbreviation chrom's name
#ifndef _FQSTATN
	static const string	Short;				// Chromosome shortening; do not convert to string in run-time
	static const chrid	UnID = -1;			// Undefined ID (255)
	static const chrid	Count = 24;			// Count of chromosomes by default (for container reserving)
	static const BYTE	MaxShortNameLength;	// Maximal length of short chrom's name
	static const BYTE	MaxNamedPosLength;	// Maximal length of named chrom's position 'chrX:12345'

private:
	static const string	sTitle;		// Chromosome title; do not convert to string in run-time
	static const char*	userChrom;	// user-defined chrom mark
	static const char*	Marks;		// heterosome marks
	static const string	UndefName;	// string not to convert in run-time

	static chrid userCID;			// user-defined chrom ID
	static chrid firstHeteroID;		// first heterosome (X,Y,M) ID
	static bool	relativeNumbering;		// true if relative numbering discipline is applied

	// Gets heterosome ID (relative or absolute) by mark,
	// or undefined ID (if absolute mode is set and cMark is inappropriate)
	static chrid HeteroID	(const char cMark);

	// Gets chromosome relative ID by mark
	static chrid GetRelativeID	(const char* cMark);

public:
	//*** common methods

	// Returns true if chromosome is autosome, false for somatic
	//	@param cad: chromosome's relative ID
	static bool IsAutosome(chrid cid) { return cid < firstHeteroID; }

	// Gets the length of prefix, or -1 if prefix is not finded
	static short PrefixLength(const char* cName);

	// Sets relative numbering discipline
	static void SetRelativeMode() { relativeNumbering = true; }

	//*** ID getters

	// Gets chromosome's ID by name
	//	@param cName: chromosome's arbitrary name
	//  @param prefixLen: length of name prefix before mark
	static chrid ID(const char* cName, size_t prefixLen=0);
	
	// Gets chrom ID by mark
	static chrid ID(const string& cMark) { return ID(cMark.c_str()); }

	// Gets chrom ID by abbreviation name
	static chrid IDbyAbbrName(const char* cAbbrName) { return ID(cAbbrName, strlen(Abbr)); }

	// Gets chrom's ID by arbitrary name
	static chrid IDbyLongName(const char* cName) { return ID(cName, PrefixLength(cName)); }

	//*** validation methods

	// Validates chromosome's arbitrary name and returns chrom ID
	//	@param cName: chromosome's arbitrary name
	//  @param prefixLen: length of prefix before chromosome's mark
	static chrid ValidateID(const char* cName, size_t prefixLen = 0);

	// Validates chromosome's mark and returns chrom ID or undefined ID
	//	@param cMark: string containing chromosome's mark
	static chrid ValidateID(const string& cMark)	{ return ValidateID(cMark.c_str()); }

	// Validates chrom abbreviation name and returns chrom ID
	//	@param cName: C-string starting with abbreviation
	static chrid ValidateIDbyAbbrName(const char* cName) { return ValidateID(cName, strlen(Abbr)); }

	// Validates all chrom ID by SAM header, applies a function to each, and sets custom ID
	//	@param samHeader: SAM header data
	//	@param f: function applied to each chromosome
	//	@param callFunc: if false do not apply function
	static void ValidateIDs(const string& samHeader, function<void(chrid cID, const char* header)> f, bool callFunc = false);

	//*** custom ID getters, setters

	// Gets ID of user specified chromosome
	static chrid UserCID()	{ return userCID; }

	// Returns true if some chromosome is specified by user
	static bool IsSetByUser() { return userCID != UnID; }
	//	@param cid: chrom id specified by user, or any chrom specified by user by default
	//static bool IsSetByUser(chrid cid = UnID) { return cid == UnID ? userCID != UnID : cid == userCID; }

	// Sets custom chromosome with control of its presence in the genome
	//	@param prColon: if true then print ": " before exception message
	//	@returns: true if custom chromosome is set successfully
	//	@throws: wrong chrom
	static bool SetUserCID(bool prColon = false);

	// Sets user-defined chromosome by mark (case-insensitive)
	static void SetUserChrom(const char* cMark);

	//*** work with name

	static const string Title(bool pl = false) { return pl ? (sTitle + 's') : sTitle;	}

	// Gets mark length by ID; used in WigMap (bioCC)
	static BYTE MarkLength(chrid cID);

	// Returns mark by ID
	static const string Mark(chrid cid);

	// Locate chromosome's mark in string.
	//	@param str: C-string checked for chromosome's abbreviation
	//	@returns: pointer to the chromosome's mark in str, or a null pointer if abbreviation is not found
	static const char* FindMark(const char* str);

	// Returns abbreviation name 'chrX'
	//	@param cid: chromosome's ID
	//	@param numbSep: if true then separate chromosome's number from the abbreviation
	static string AbbrName(chrid cid, bool numbSep = false);

	// Returns short name 'chrom X'
	//	@param cid: chromosome's ID
	static string ShortName(chrid cid)	{ return Short + SPACE + Mark(cid); }

	// Returns title name 'chromosome X' or 'chromosomes'
	//	@param cid: chromosome's ID or UnID if plural
	static const string TitleName(chrid cid = UnID);

	// Returns 'there is no chromosome <user-defined-chrom>' string
	static const string NoChromMsg();

	//static const string Absent(chrid cid, const string& what);

#endif	// _FQSTATN
} chrom;


// 'Region' represents a simple region within nucleotides array (chromosome).
struct Region
{
	chrlen Start;	// start position of the region in standard chromosomal coordinates
	chrlen End;		// end position of the region in standard chromosomal coordinates

	Region(chrlen start = 0, chrlen end = 0) : Start(start), End(end) {}

	// Constructs Read
	//	@param r: original read
	Region(const Region& r) { memcpy(this, &r, sizeof(Region)); }

	// Constructs extended Read (fragment)
	//	@param r: original read
	//	@param extLen: extended length
	//	@param reverse: true if read is reversed (neg strand)
	Region(const Region& r, fraglen extLen, bool reverse);

	// Initializes instance
	void Set(chrlen start = 0, chrlen end = 0) { Start = start; End = end; }

	chrlen Length()	const { return End - Start; } // The End is not included in the bases https://genome.ucsc.edu/FAQ/FAQformat.html#format1

	bool Empty()	const { return !End; }

	chrlen Centre()	const { return Start + (Length() >> 1); }

	bool operator==(const Region& r) const { return !(*this != r); }
	bool operator!=(const Region& r) const { return memcmp(this, &r, sizeof(Region)); }

	void operator-=(chrlen val) { Start -= val; End -= val; }

	// Returns true if this instance is invalid
	bool Invalid() const { return Start >= End; }

	// Returns true if Region r is covered by this instance.
	//	@param r: tested Region; should be sorted by start position
	bool PlainCover(const Region& r) const { return r.End <= End && r.Start >= Start; }

	// Returns true if Region r is adjoined with this instance.
	//	@param r: tested Region; should be sorted by start position
	bool Adjoin(const Region& r) const { return r.Start == End; }

	// Returns true if Region r is crossed with this instance.
	//	@param r: tested Region; should be sorted by start position
	bool Cross(const Region& r) const { return r.Start < End && r.End > Start; }

	// Compares two Regions by start position. For sorting a container.
	static bool CompareByStartPos(const Region& r1, const Region& r2) {	return r1.Start < r2.Start;	}

	// Extends Region with chrom length control.
	// If extended Region starts from negative, or ends after chrom length, it is fitted.
	//	@param extLen: extension length in both directions
	//	@param cLen: chrom length; if 0 then no check
	void Extend(chrlen extLen, chrlen cLen);

	void Print() const { cout << Start << TAB << End << LF; }

private:
	void ExtStart(fraglen extLen) { Start = End - extLen; }
	void ExtEnd	 (fraglen extLen) { End = Start + extLen; }

	typedef void (Region::* tExtStartEnd)(fraglen);
	static tExtStartEnd fExtStartEnd[2];
};

// 'Regions' represents a container of defined regions within chromosome
// Defined in TxtFile.h since it's used in Fa class
class Regions
{
protected:
	vector<Region> _regions;

public:
	typedef vector<Region>::const_iterator Iter;

	// Iterator to region Begin position
	const Iter Begin()	const { return _regions.begin(); }

	// Iterator to region End position
	const Iter End()		const { return _regions.end(); }

	// Default (empty) constructor to form Chroms collection
	Regions() {}

	// Single region constructor
	Regions(chrlen start, chrlen end) { _regions.emplace_back(start, end); }

	// Copying constructor
	//Regions(const Regions& rgns) { _regions = rgns._regions;	}

	// Gets total length of regions.
	//chrlen Length() const;

	// Gets count of regions.
	chrlen Count()		const { return chrlen(_regions.size()); }

	// Gets first start position.
	chrlen FirstStart()	const { return _regions.front().Start; }

	// Gets last end position.
	chrlen LastEnd()	const { return _regions.back().End; }

	// Gets conditionally defined length: distance between first start and last end
	chrlen DefLength()	const { return LastEnd() - FirstStart(); }

	//Regions& operator=(const Regions& rgn);

	const Region& operator[](chrlen ind) const { return _regions[ind]; }

	// Reserves container's capacity.
	//	@capacity: reserved number of regions. The real number may be differ.
	void Reserve(size_t capacity) { _regions.reserve(capacity); }

	// Clears container.
	void Clear() { _regions.clear(); }

	void Add(const Region& rgn) { _regions.push_back(rgn); }

	void Add(chrlen start, chrlen end) { _regions.emplace_back(start, end); }

	// Copies subregions
	void Copy(const vector<Region>& source, chrlen start, chrlen stop) {
		_regions = vector<Region>(source.begin() + start, source.begin() + stop + 1);
	}

#if defined _READDENS || defined _BIOCC

	// Returns an iterator referring to the past-the-end element, where end is external
	//	@curr_it: region's const iterator, from which the search is started
	//	@end: external pre-defined end coordinate
	Iter ExtEnd(Iter curr_it, chrlen end) const;

	// Initializes this instance by intersection of two Regions.
	void FillOverlap(const Regions& regn1, const Regions& regn2);

	// Initializes this instance by inverted external Regions.
	//	@regn: external Regions
	//	@masEnd: the maximum possible end-coordinate of region:
	//	the chromosome length in case of nucleotides sequance.
	void FillInvert(const Regions& regn, chrlen maxEnd);

	// Initializes this instance by external Regions.
	void Copy(const Regions& rgns) { _regions = rgns._regions; }

#endif	// _READDENS, _BIOCC
#ifdef DEBUG
	void Print() const;
#endif
};

#endif	// _COMMON_H