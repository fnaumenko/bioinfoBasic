/**********************************************************
common.cpp
Last modified: 04/30/2024
***********************************************************/

#include "common.h"
#include <sstream>
#include <algorithm>	// std::transform, REAL_SLASH
#ifdef OS_Windows
#define SLASH '\\'		// standard Windows path separator
#define REAL_SLASH '/'	// is permitted in Windows too
#else
#define SLASH '/'		// standard Linux path separator
#endif

/************************ common Functions ************************/

chrlen atoui(const char* p)
{
	chrlen x = 0;
	for (; *p >= '0' && *p <= '9'; ++p)
		x = x * 10 + (*p - '0');
	return x;
}

size_t atoul(const char* p)
{
	size_t x = 0;
	for (; *p >= '0' && *p <= '9'; ++p)
		x = x * 10 + (*p - '0');
	return x;
}

int OnesCount(int n)
{
	int cnt = 0;
	//for (; n; n >>= 1)  cnt += n & 1;		// sample 11001: 5 cycles
	for (; n; n &= n - 1)   cnt++;			// sample 11001: 3 cycles
	return cnt;
}

int RightOnePos(int n)
{
	int pos = 0;
	for (; n ^ 1; n >>= 1)	pos++;
	return pos;
}

int DigitsCount(size_t val, bool isLocale)
{
	int res = 0;
	for (; val; val /= 10, res++);
	if (isLocale)	res += (res - 1) / 3;
	return res;
}

string	PercentToStr(float val, BYTE precision, BYTE fieldWith, bool parentheses)
{
	//const float minVal = float(pow(10.f, -precision));	// minimum printed value
	const float minVal = 1.f / float(pow(10.f, precision));	// minimum printed value
	stringstream ss;

	ss << SPACE;
	if (parentheses)		ss << '(';
	//if(fieldWith)	ss << setw(--fieldWith);	// decrease to account '%' sign
	if (val && val < minVal) {
		if (fieldWith) {
			int blankCnt = fieldWith - precision - 5;
			if (blankCnt > 0)	ss << setw(blankCnt) << SPACE;
		}
		ss << '<' << minVal;
	}
	else {
		if (precision && val) {
			if (val >= 100)	precision = 3;
			else if (val < 1)	ss << fixed;
			ss << setprecision(precision);
		}
		//if(val - int(val) != 0)	
		//	ss << fixed << setprecision(precision);
		if (fieldWith)	ss << setw(fieldWith -= 2);	// decrease one for SPACE, one for '%'
		ss << val;
	}
	ss << PERS;
	if (parentheses)		ss << ')';
	return ss.str();
}

void PrintHorLine(int lw)
{
#ifdef OS_Windows
	wcout << setw(++lw) << setfill(L'\304') << L'\n';
#ifdef _DUP_OUTPUT
	dout.File() << string(lw, HPH) << LF;
#endif
#else
	for (int i = 0; i < lw; dout << "─", i++);	dout << LF;
	//dout << setw(++lw) << setfill('─') << LF;
#endif
}

#if defined _WIGREG || defined _BIOCC

chrlen AlignPos(chrlen pos, BYTE res, BYTE relative)
{
	short rest = pos % res;
	return rest - relative ?
		pos + relative - rest + (rest > res << 1 ? res : 0) :
		pos;
}

#endif

//#ifdef OS_Windows
//string Wchar_tToString(const wchar_t* wchar)
//{
//	string str = strEmpty;
//	while(*wchar)	str += char(*wchar++);
//	return str;
//}
//
//wchar_t* StringToWchar_t(const string &str, wchar_t* wchar)
//{
//	int i = (int)str.size();
//	wchar[i] = 0;
//	for(--i; i >= 0; i--)
//		wchar[i] = (wchar_t)str[i];
//	return wchar;
//}
//#endif

/************************ end of common Functions ************************/

#ifdef _DUP_OUTPUT

bool dostream::OpenFile(const string fname)
{
	if (!fname.length())	return false;
	file.open(fname.c_str());
	return true;
}

#endif

PairVals::PairVals(float val1, float val2, float min1, float min2, float max1, float max2)
{
	vals[SET] = make_pair(val1, val2);
	vals[MIN] = make_pair(min1, min2);
	vals[MAX] = make_pair(max1, max2);
}

/************************ class Err ************************/

//const char* Err::TREAT_BED_EXT = "after extension";	// clarifying message in the bed stretch operation

const char* Err::_msgs[] = {
	/* NONE */		"WARNING",
	/* MISSED */	"missing",
	/* F_NONE */	"no such file",
	/* D_NONE */	"no such folder",
	/* FD_NONE */	"no such file or folder",
	/* F_MEM */		"memory exceeded",
	/* F_OPEN */	"could not open",
	/* F_CLOSE */	"could not close",
	/* F_READ */	"could not read",
	/* F_EMPTY */	"empty",
	/* F_BIGLINE */	"record length exceeds buffer limit",
	/* FZ_MEM */	"not enough internal gzip buffer",
	/* FZ_OPEN */	"wrong reading mode READ_ANY for gzip file",
	/* FZ_BUILD */	"this build does not support gzip files",
	/* F_WRITE */	"could not write",
	///* F_FORMAT */	"wrong format",
	#ifndef _FQSTATN
	/* TF_FIELD */	"number of fields is less than required",
	/* TF_EMPTY */	"no records",
	#endif
	/* EMPTY */		""
};

const char* Err::FailOpenOFile = "could not open output file";

void StrCat(char* dst, const char* src) { memcpy(dst + strlen(dst), src, strlen(src) + 1); }

void Err::set_message(const char* sender, const char* txt, const char* specTxt)
{
	size_t senderLen = sender != NULL ? strlen(sender) : 0;
	size_t outLen = senderLen + strlen(txt) + strlen(SepSCl) + 2;
	if (specTxt)	outLen += strlen(specTxt) + 1;
	_outText = new char[outLen];		// will be free in destructor
	memset(_outText, cNULL, outLen);
	if (sender) {
		//*_outText = SPACE;
		if (senderLen)	StrCat(_outText, sender);
		StrCat(_outText, SepCl);
	}
	StrCat(_outText, txt);
	if (specTxt) {
		if (*specTxt != ':')		StrCat(_outText, sSPACE);
		StrCat(_outText, specTxt);
	}
}

//const string Err::IssueNumbToStr(const string& issName, ULONG issNumb, const string& fName)
//{
//	string res = fName;
//	if(fName != strEmpty)	res += SepSCl;
//	return res + issName + SPACE + NSTR(issNumb);
//}

const string Err::MsgNoFile(const string& fName, bool plural, const string fExt)
{
	const string pl_ext = plural ? "s" : strEmpty;
	return string("no " + fName + fExt + "[" + ZipFileExt + "] file" + pl_ext + " in this directory");
}


Err::Err(const Err& src)
{
	_code = src._code;
	auto size = strlen(src.what()) + 1;
	_outText = new char[size];
	memset(_outText, 0, size);
	memcpy(_outText, src._outText, size - 1);	// instead of strcpy(_outText, src._outText)
}

void Err::Throw(bool throwExc, bool eol) {
	if (throwExc)
		throw* this;
	else {
		dout << what();
		if (eol)		dout << LF;
		fflush(stdout);
	}
}

void Err::Warning(bool eol, bool prefix)
{
	if (prefix)	dout << SepCl;
	else		dout << LF;
	dout << _msgs[ErrWARNING];
	if (*what() != ':') dout << SepCl;	// check if sender is not recorded
	dout << what();
	if (eol)	dout << LF;
	fflush(stdout);
}

/************************ end of class Err ************************/

/************************ class FileSystem ************************/

bool FS::IsExist(const char* name, int st_mode)
{
	struct_stat64 st;
	const auto len = strlen(name) - 1;
	int res;

	if (name[len] == SLASH) {
		string sname = string(name, len);
		res = _stat64(sname.c_str(), &st);
	}
	else res = _stat64(name, &st);
	return (!res && st.st_mode & st_mode);
	//return ( !_stat64(name, &st) && st.st_mode & st_mode );
}

bool FS::CheckExist(const char* name, int st_mode, bool throwExcept, Err::eCode ecode)
{
	if (IsExist(name, st_mode))	return false;
	Err(ecode, name).Throw(throwExcept);
	return true;
}

size_t FS::GetLastExtPos(const string& fname) {
	size_t pos = fname.find_last_of(DOT);
	return (pos != string::npos && (pos == 0 || (pos == 1 && fname[0] == DOT))) ?	// ./name || ../name
		string::npos : pos;
}

bool SearchExt(const string& fname, const string& ext, bool isZip, bool composite)
{
	size_t pos = fname.find(ext);
	if (pos == string::npos)		return false;
	if (!composite)				return fname.size() - pos == ext.size();
	return fname.size() - pos - isZip * ZipFileExt.length() == ext.size();
}

bool FS::HasCaseInsExt(const string& fname, const string& ext, bool knownZip, bool composite)
{
	bool res = SearchExt(fname, ext, knownZip, composite);
	if (!res) {		// try to case insensitive search
		string str(fname);
		string substr(ext);
		transform(str.begin(), str.end(), str.begin(), ::tolower);
		transform(substr.begin(), substr.end(), substr.begin(), ::tolower);
		res = SearchExt(str, substr, knownZip, composite);
	}
	return res;
}


LLONG FS::Size(const char* fname)
{
	struct_stat64 st;
	return _stat64(fname, &st) == -1 ? -1 : st.st_size;
}
//// alternative implementation in C++ style;
// needs #include <fstream>
//std::ifstream::pos_type filesize(const char* fname)
//{
//    std::ifstream ifs(fname, std::ifstream::ate | std::ifstream::binary);
//    return ifs.tellg(); 
//}

LLONG FS::UncomressSize(const char* fname)
{
	FILE* file = fopen(fname, "rb");	// "read+binary"
	if (file == NULL)		return -1;
	BYTE sz[] = { 0,0,0,0 };
	_fseeki64(file, -4, SEEK_END);
	fread(sz, 1, 4, file);
	fclose(file);
	return (sz[3] << 3 * 8) + (sz[2] << 2 * 8) + (sz[1] << 8) + sz[0];
}

bool FS::CheckFileDirExist(const char* name, const string& ext, bool throwExcept)
{
	return HasExt(name, ext) ?
		CheckFileExist(name, throwExcept) :
		CheckFileDirExist(name, throwExcept);
}

const char* FS::CheckedFileDirName(const char* name)
{
	if (!IsFileDirExist(name))	Err(Err::FD_NONE, name).Throw();
	return name;
}

const char* FS::CheckedFileName(const char* name)
{
	if (name && !IsFileExist(name))	Err(Err::F_NONE, name).Throw();
	return name;
}

//const char* FS::CheckedDirName(int opt)
//{
//	const char* name = Options::GetSVal(opt);
//	if (name && !IsDirExist(name))	Err(Err::D_NONE, name).Throw();
//	return name;
//}

string const FS::GetExt(const char* fname) {
	const char* pdot = strrchr(fname, DOT);
	if (!pdot)		return strEmpty;
	if (strcmp(pdot, ZipFileExt.c_str()))
		return string(pdot + 1);				// no zip extention
	const char* pprevdot = pdot - 1;
	for (; pprevdot >= fname; pprevdot--)	// find previous DOT
		if (*pprevdot == DOT)
			break;
	return pprevdot + 1 == fname ? strEmpty : string(pprevdot + 1, pdot - pprevdot - 1);
}

bool FS::IsShortFileName(const string& fname)
{
#ifdef OS_Windows
	if (fname.find_first_of(REAL_SLASH) != string::npos)	return false;
	else
#endif
		return fname.find_first_of(SLASH) == string::npos;
}

string const FS::ShortFileName(const string& fname)
{
#ifdef OS_Windows
	if (fname.find(REAL_SLASH) != string::npos) {
		string tmp(fname);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp.substr(tmp.find_last_of(SLASH) + 1);
	}
#endif
	return fname.substr(fname.find_last_of(SLASH) + 1);
}

string const FS::DirName(const string& fname, bool addSlash)
{
#ifdef OS_Windows
	if (fname.find(REAL_SLASH) != string::npos) {
		string tmp(fname);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp.substr(0, tmp.find_last_of(SLASH) + int(addSlash));
	}
#endif
	return fname.substr(0, fname.find_last_of(SLASH) + int(addSlash));
}

string const FS::LastSubDirName(const string& fname)
{
	const string& dir = FS::DirName(fname, false);
	size_t pos = dir.find_last_of(SLASH);
	return pos == string::npos ? dir :
		dir.substr(dir.substr(0, pos).length() + 1);
}

string const FS::LastDirName(const string& name)
{
	size_t pos = name.find_last_of(SLASH);
	if (pos == string::npos)		return name;
	if (++pos == name.length()) {				// ended by slash?
		pos = name.find_last_of(SLASH, pos - 2) + 1;
		return name.substr(pos, name.length() - 1 - pos);
	}
	return name.substr(pos);
}

string const FS::MakePath(const string& name)
{
#ifdef OS_Windows
	if (name.find(REAL_SLASH) != string::npos) {
		string tmp(name + SLASH);
		replace(tmp.begin(), tmp.end(), REAL_SLASH, SLASH);
		return tmp;
	}
#endif
	return name[name.length() - 1] == SLASH ? name : name + SLASH;
}

#if !defined _WIGREG && !defined _FQSTATN
bool FS::GetFiles(vector<string>& files, const string& dirName, const string& ext, bool all)
{
#ifdef OS_Windows
	int count = 0;
	string fileTempl = FS::MakePath(dirName) + '*' + ext;
	WIN32_FIND_DATAA ffd;

	HANDLE hFind = FindFirstFileA(LPCSTR(fileTempl.c_str()), LPWIN32_FIND_DATAA(&ffd));
	if (hFind == INVALID_HANDLE_VALUE)
		return false;
	if (all) {
		// count files to reserve files capacity
		do	count++;
		while (FindNextFileA(hFind, LPWIN32_FIND_DATAA(&ffd)));
		files.reserve(count);
		// fill files
		hFind = FindFirstFileA(LPCSTR(fileTempl.c_str()), LPWIN32_FIND_DATAA(&ffd));
		do	files.emplace_back(ffd.cFileName);
		while (FindNextFileA(hFind, LPWIN32_FIND_DATAA(&ffd)));
	}
	else
		files.emplace_back(ffd.cFileName);
	FindClose(hFind);
	return true;
#else
	struct dirent* entry;
	DIR* dir = opendir(dirName.c_str());	// doesn't need to check because of appl options inspection
	string name;

	// fill files
	files.reserve(Chrom::Count);
	while (entry = readdir(dir)) {
		name = string(entry->d_name);
		if (HasExt(name, ext, false)) {
			files.push_back(name);
			if (!all)	break;
		}
	}
	closedir(dir);
	return files.size() > 0;
#endif	// OS_Windows
}
#endif	// _WIGREG, _FQSTATN
/************************ end of class FileSystem ************************/

/************************  class TimerBasic ************************/

bool	TimerBasic::Enabled = false;

void PrintTime(long elapsed, bool parentheses, bool isLF)
{
	using namespace std::chrono;

	//if (offset)	cout << setw(offset) << setfill(' ') << ' ';
	const auto hrs = uint16_t(elapsed / 1000 / 60 / 60);
	const auto mins = uint16_t(elapsed / 1000 / 60 % 60);
	const auto secs = uint16_t(elapsed / 1000 % 60);
	const bool prMins = mins || secs > 9;
	auto Round = [](int val, int iterations) {
		int result = val;
		for (auto i = 0; i < iterations; i++) {
			int rem = result % 10;
			result -= rem;
			if (rem >= 5)	result += 10;
			result /= 10;
		}
		return result;
	};

	if (parentheses)	dout << '(';
	dout << setfill('0') << right;		// right couse it may be chanched by previuos output
	//cout.width(2);	applied to the first 'cout' only...
	if (hrs)	dout << setw(2) << hrs << COLON;
	if (prMins)	dout << setw(2) << mins << COLON;
	dout << setw(2) << secs;
	if (!hrs && !prMins)
		dout << '.' << Round(elapsed % 1000, 1 + bool(secs));
	if (parentheses)	dout << ')';
	if (isLF)	dout << LF, fflush(stdout);
}

long TimerBasic::GetElapsed() const
{
	using namespace std::chrono;

	return long(duration_cast<milliseconds>(steady_clock::now() - _begin).count());
}

void TimerBasic::Print(long elapsed, const char* title, bool parentheses, bool isLF)
{
	if (title)	dout << title;
	PrintTime(elapsed, parentheses, isLF);
}

/************************  end ofclass TimerBasic ************************/

/************************  class Timer ************************/
clock_t	Timer::_StartCPUClock;

void Timer::Stop(int offset, bool parentheses, bool isLF)
{
	if (_enabled) {
		if (offset)	dout << setw(offset) << SPACE;
		PrintTime(GetElapsed(), parentheses, isLF);
	}
}

/************************  end of class Timer ************************/
#ifdef _TEST
/************************  class Stopwatch ************************/

void Stopwatch::Stop(const string title) const
{
	if (_isStarted) {
		_sumTime += GetElapsed();
		if (title != strEmpty)	PrintTime(_sumTime, false, true);
	}
}

/************************  end of class Stopwatch ************************/

/************************  class StopwatchCPU ************************/

void StopwatchCPU::Stop(const char* title, bool print, bool isLF)
{
	_sumclock += clock() - _clock;
	if (print) {
		if (title)	dout << title;
		PrintTime(1000 * _sumclock / CLOCKS_PER_SEC, false, isLF);
	}
}

/************************  end of class StopwatchCPU ************************/
#endif	// _TEST
#ifdef _MULTITHREAD
/************************  class Mutex ************************/

bool	Mutex::_active;		// true if the mutex really should work
mutex	Mutex::_mutexes[int(Mutex::eType::NONE)];

/************************  end of class Mutex ************************/
#endif	// _MULTITHREAD

/************************ class Chrom ************************/

const char*		Chrom::Abbr = "chr";
const BYTE		Chrom::MaxAbbrNameLength = 5;
#ifndef _FQSTATN
const char cM = 'M';
const char*		Chrom::userChrom;
const char*		Chrom::Marks = "XYM";
const string	Chrom::UndefName = "UNDEF";
const string	Chrom::Short = "chrom";
const string	Chrom::sTitle = "chromosome";
const BYTE		Chrom::MaxShortNameLength = BYTE(Short.length()) + MaxMarkLength;
const BYTE		Chrom::MaxNamedPosLength = BYTE(strlen(Abbr)) + MaxMarkLength + CHRLEN_CAPAC + 1;


chrid Chrom::userCID = UnID;
chrid Chrom::firstHeteroID;
bool Chrom::relativeNumbering;

chrid Chrom::HeteroID(const char cMark)
{
	if (!relativeNumbering)
		// absolute heterosome ID; in case of 'M' multiply by 2 so as not to violate the sorting check in input BAM/BED
		return cMark << int(cMark == cM);
	for (size_t i = 0; i < strlen(Marks); i++)
		if (cMark == Marks[i])
			return firstHeteroID + chrid(i);	// relative heterosome ID
	return UnID;
}

chrid Chrom::GetRelativeID(const char* cMark)
{
	if (isdigit(*cMark)) {		// autosome
		//chrid id = atoi(cMark) - 1;
		chrid id = chrid(atoui(cMark) - 1);
		return id < firstHeteroID ? id : UnID;
	}
	return HeteroID(*cMark);	// heterosome
}

const char* SubStr(const char* str, const char* templ, size_t templLen)
{
	str = strchr(str, *templ);
	if (str)
		for (size_t i = 1; i < templLen; i++)
			if (*++str != templ[i])
				return SubStr(str, templ, templLen);
	return str;
}

short Chrom::PrefixLength(const char* cName)
{
	// search from the beginning of the cName because
	// it simpler to find the first digit for multidigit mark

	// start search from the occurrence of 'chr'
	const char* str = SubStr(cName, Abbr, strlen(Abbr));
	if (str)
		for (; *str; str++)
			if (isdigit(*str) || isupper(*str))
				return short(str - cName);
	return -1;
}

chrid Chrom::ID(const char* cName, size_t prefixLen)
{
	//return isdigit(*(cName += prefixLen)) ? atoi(cName) - 1 : HeteroID(*cName);
	return isdigit(*(cName += prefixLen)) ? chrid(atoui(cName)) - 1 : HeteroID(*cName);
}

chrid Chrom::ValidateID(const char* cName, size_t prefixLen)
{
	if (!cName)					return UnID;
	cName += prefixLen;							// skip prefix
	for (int i = 1; i <= MaxMarkLength; i++)
		if (cName[i] == USCORE)	return UnID;	// exclude chroms with '_'

	if (isdigit(*cName)) {						// autosome
		//chrid id = atoi(cName);
		chrid id = atoui(cName);
		if (relativeNumbering && id > firstHeteroID)	firstHeteroID = id;
		return id - 1;
	}
	return HeteroID(*cName);					// heterosome
}

void Chrom::ValidateIDs(const string& samHeader, function<void(chrid cID, const char* header)> f, bool callFunc)
{
	relativeNumbering = true;
	for (const char* header = samHeader.c_str();
		header = strstr(header, Abbr);
		header = strchr(header, LF) + strlen("\n@SQ\tSN:"))
	{
		chrid cID = ValidateIDbyAbbrName(header);
		if (callFunc)
			f(cID, strchr(header, TAB) + strlen("\tLN:"));
	}
	SetUserCID(true);
}

void Chrom::SetUserCID(bool prColon)
{
	if (userChrom && (userCID = GetRelativeID(userChrom)) == UnID)
		Err((prColon ? SepCl : strEmpty) + NoChromMsg() + " in this genome").Throw();
}

void Chrom::SetUserChrom(const char* cMark)
{
	if (cMark) {
		userChrom = cMark;
		*(const_cast<char*>(userChrom)) = toupper(*cMark);
		userCID = ValidateID(userChrom);
	}
}

BYTE Chrom::MarkLength(chrid cID)
{
	return cID == UnID ? BYTE(UndefName.length()) : (cID >= firstHeteroID || cID < 9 ? 1 : 2);
}

const string Chrom::Mark(chrid cid)
{
	auto autosomeToStr = [](chrid cid) { return to_string(cid + 1); };

	if (cid == UnID)		return UndefName;
	if (relativeNumbering)
		return cid < firstHeteroID ? autosomeToStr(cid) :
		(cid > firstHeteroID + 2) ? UndefName : string{ Marks[cid - firstHeteroID] };
	return cid < '9' ? autosomeToStr(cid) : string{ char(cid >> int(cid == 2*cM)) };	// divide by 2 in case 'M'
}

const char* Chrom::FindMark(const char* str)
{
	const char* substr = strstr(str, Abbr);
	return substr ? (substr + strlen(Abbr)) : NULL;
}

string Chrom::AbbrName(chrid cid, bool numbSep)
{
	return Abbr + (numbSep ? sSPACE : strEmpty) + Mark(cid);
}

const string Chrom::TitleName(chrid cid)
{
	return sTitle + (cid == UnID ? "s" : SPACE + Mark(cid));
}

const string Chrom::NoChromMsg()
{
	// Checking for 'userChrom' is redundant, because checking for empty sequence is already done in the 'UniBedReader' constructor.
	// It was made for testing purposes and doesn't affect performance, because this is the final output
	return "there is no " + (userChrom ? sTitle +sSPACE + userChrom : Title(true));
}

//const string Chrom::Absent(chrid cid, const string& what)
//{
//	return AbbrName(cid) + " is absent in " + what + " file: skipped";
//}

#endif	// _FQSTATN

/************************ end of class Chrom ************************/

/************************ struct Region ************************/

Region::tExtStartEnd Region::fExtStartEnd[2] = { &Region::ExtEnd, &Region::ExtStart };

Region::Region(const Region& r, fraglen extLen, bool reverse)
{
	memcpy(this, &r, sizeof(Region));
	(this->*fExtStartEnd[reverse])(extLen);
}

void Region::Extend(chrlen extLen, chrlen cLen)
{
	Start -= extLen > Start ? Start : extLen;
	End += extLen;
	if (cLen && End > cLen)	End = cLen;
}


/************************ class Regions ************************/

//chrlen Regions::Length () const
//{
//	chrlen len = 0;
//	for(Iter it=_regions.begin(); it<_regions.end(); it++)
//		len += it->Length();
//	return len;
//}

//Regions& Regions::operator=(const Regions& rgn)
//{
//	_regions = rgn._regions;
//	return *this;
//}

#if defined _READDENS || defined _BIOCC

Regions::Iter Regions::ExtEnd(Iter curr_it, chrlen end) const
{
	Iter it = curr_it;
	for (; it != _regions.end(); it++)
		if (it->End > end)	break;
	return it;
}

void Regions::FillOverlap(const Regions& regn1, const Regions& regn2)
{
	chrlen start = 0, end = 0, start1, start2, end1, end2;
	Iter it1 = regn1._regions.begin();
	Iter it2 = regn2._regions.begin();
	Reserve(max(regn1.Count(), regn2.Count()));
	for (; it1 != regn1._regions.end(); it1++) {
		start1 = it1->Start;	end1 = it1->End;
		start2 = it2->Start;	end2 = it2->End;
		if (start1 < start2) {
			if (end1 > start2) {
				start = start2;
				if (end1 > end2) { end = end2;	it2++; it1--; }
				else				end = end1;
			}
		}
		else
			if (start1 >= end2) { it2++; it1--; }
			else {
				start = max(start1, start2);
				if (end1 > end2) { end = end2;	it2++; it1--; }
				else				end = end1;
			}
		if (end) {
			Add(start, end);
			if (it2 == regn2._regions.end())		return;
			end = 0;
		}
	}
}

void Regions::FillInvert(const Regions& regn, chrlen maxEnd)
{
	Region rgn;
	Iter it = regn._regions.begin();

	Reserve(regn.Count() + 1);
	for (; it != regn._regions.end(); it++) {
		rgn.End = it->Start - 1;
		Add(rgn);
		rgn.Start = it->End + 1;
	}
	if (rgn.Start < (maxEnd + 1)) {
		rgn.End = maxEnd;
		Add(rgn);
	}
}

#endif	// _READDENS, _BIOCC

#ifdef DEBUG
void Regions::Print() const
{
	int i = 0;
	cout << "Regions:\n";
	for (Iter it = _regions.begin(); it < _regions.end(); it++)
		cout <<
		//setw(2) << 
		++i << COLON << TAB <<
		//setw(9) << 
		it->Start << TAB << it->End << endl;
}
#endif
/************************ end of class Regions ************************/
