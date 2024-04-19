/**********************************************************
common.cpp
Last modified: 04/18/2024
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
	dout.ToFile(string(lw, HPH) + LF);
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

/************************ class Options ************************/

#define ENUM_REPLACE '?'	// symbol in description that is replaced by enum value

#define PRINT_IN_PRTHS(v)	cout<<" ["<<(v)<<']'	// prints value in parentheses
#define OPT_TO_STREAM(opt)	HPH<<(opt)

const char* optTitle = "option ";
const char* Spotteruous = "Spotteruous";
const char* Default = " Default: ";
const char* Missing = "missing ";
const char* Warning = "WARNING: ";
const string sValue = "value";

const char* Options::sPrSummary = "print program's summary";
const char* Options::sPrTime = "print run time";
const char* Options::sPrUsage = "print usage information";
const char* Options::sPrVersion = "print program's version";
const char* Options::Booleans[] = { "OFF","ON" };
const char* Options::TypeNames[] = {
	NULL, "<name>", "<char>", "<int>", "<float>", "<long>", NULL, NULL,
	"<[int]:[int]>", "<[float]:[float]>", NULL, NULL, NULL
};

const char Options::Option::EnumDelims[] = { '|', ',', ':' };

bool Options::Option::IsValidFloat(const char* str, bool isInt, bool isPair)
{
	char c = *str;
	const char* str0 = str;
	BYTE dotCnt = 0, eCnt = 0;

	if (!isdigit(c))						// check first char
		if (c == DOT)	dotCnt++;
		else if (c != HPH && c != PLUS)	goto end;
	for (str++; *str; str++)				// check next chars
		if ((c = *str) == DOT)
			if (dotCnt)	goto end;		// more than one dot
			else		dotCnt++;
		else if (tolower(c) == 'e')
			if (eCnt)	goto end;		// more than one 'e'
			else		eCnt++;
		else if (!isdigit(c))
			if (isPair && c == EnumDelims[2])	break;	// don't check the substring after ';'
			else goto end;				// wrong symbol
	if (isInt && dotCnt && !eCnt)		// e.g. 1.5e1 is acceptable as int
		cerr << Warning << ToStr(false) << SepSCl << "float value "
		<< (isPair ? "in " : strEmpty)
		<< str0 << " will be treated is integer\n";
	return true;
end:
	return !PrintWrong(str0);
}

void Options::Option::PrintTransformDescr(char* buff, const char** vals, short* cnt)
{
	if (vals) {			// enum/combi?
		const char* subStr = strchr(buff, ENUM_REPLACE);
		if (subStr) {	// something to replace by enum value
			size_t strLen = subStr - buff;
			buff[strLen] = 0;
			cout << buff << vals[(*cnt)++];	// output substring and enum value
			buff[strLen] = ENUM_REPLACE;
			PrintTransformDescr(buff + strLen + 1, vals, cnt);
		}
		else	cout << buff;
	}
	else	cout << buff;
}

void Options::Option::PrintSubLine(char* buff, const char* str, const char* subStr, const char** vals, short* cnt)
{
	if (subStr) {	// is substring ended by LF exist?
		// form substring
		size_t strLen = subStr - str;
		memcpy(buff, str, strLen);	// instead of strncpy(buff, str, strLen);
		buff[strLen] = 0;
		PrintTransformDescr(buff, vals, cnt);	// output enum values
		cout << LF;
		for (BYTE t = 0; t < IndentInTabs; t++)	cout << TAB;
		str = subStr + 1;		// skip LF
		subStr = strchr(str, LF);
		PrintSubLine(buff, str, subStr, vals, cnt);
	}
	else {			// output rest of initial string without LF
		memcpy(buff, str, strlen(str)); // instead of strcpy(buff, str);
		PrintTransformDescr(buff, vals, cnt);
	}
}

int Options::Option::SetVal(const char* opt, bool isword, char* val, char* nextItem, int& argInd)
{
	if (isword) {
		if (!Str || strcmp(Str, opt))	return -1;	// not this option
		Sign.Set(tOpt::WORD);
	}
	else if (Char != *opt)	return -1;		// not this option

	if (Sign.Is(tOpt::TRIMMED))	return PrintAmbigOpt(opt, isword, "duplicated");

	//== check actual value
	const bool isValOblig = ValRequired() && !IsValEsc();
	const bool noRealVal =		// true if no actual value is submitted
		val == NULL								// no value
		|| (!nextItem && !isValOblig)			// nextItem is apps parameter
		|| (*val == '-' && !isdigit(val[1]));		// val is not negative value but next option

	if (noRealVal) {
		if (isValOblig)	return PrintWrong(NULL, sValue + " required");
		val = NULL;				// in case of 'false' value
	}
	else
		if (ValRequired())	argInd++;
		else
			if (val && nextItem && *nextItem == HPH) {	// next token after val is option
				cerr << Warning, PrintWrong(NULL, sValue + " prohibited: " + string(val) + " ignored");
				argInd++;
			}

	//== set actual value
	Sign.Set(tOpt::TRIMMED);
	switch (ValType) {
	case tNAME: SVal = noRealVal ? NULL : val;	return 0;
	case tENUM:	return SetEnum(val);
	case tCOMB:	return SetComb(val);
	case tCHAR:	if (NVal != NO_VAL && strlen(val) == 1)
		return SetTriedFloat(*val, MinNVal, MaxNVal);	// value is treated as int
		break;
	case tINT:
	case tFLOAT:
	case tLONG:	return noRealVal ?
		0 :				// NVal is default val
		IsValidFloat(val, ValType == tINT) ?
		SetTriedFloat(float(atof(val)), MinNVal, MaxNVal) :
		1;
	case tHELP:	return PrintUsage(true);
	case tSUMM:	return PrintSummary(false);
	case tVERS:	return PrintVersion();
	default:	return SetPair(noRealVal ? NULL : val, ValType == tPR_INT);	// tPR_INT, tPR_FL
	}
	return PrintWrong(val);
}

int Options::Option::CheckOblig() const
{
	if (Sign.Is(tOpt::OBLIG) && ValRequired()
		&& ((ValType == tNAME && !SVal) || (ValType != tNAME && NVal == NO_VAL))) {
		cerr << Missing << "required option " << NameToStr(false) << LF;
		return -1;
	}
	return 1;
}

string Options::Option::NameToStr(bool asPointed) const
{
	ostringstream ss;
	ss << HPH;
	if (asPointed) {
		if (!Sign.Is(tOpt::WORD)) { ss << Char; return ss.str(); }
	}
	else if (Char != HPH) {
		ss << Char;
		if (Str)	ss << '|' << HPH;
	}
	if (Str)		ss << HPH << Str;
	return ss.str();
}

const string Options::Option::PairValsToStr(const pairVal* vals) const
{
	static const char* sAuto = "auto";
	ostringstream ss;
	ss << dec;
	if (vals->first == vUNDEF)	ss << sAuto;
	else						ss << vals->first;
	ss << EnumDelims[2];
	if (vals->second == vUNDEF)	ss << sAuto;
	else						ss << vals->second;
	return ss.str();
	//return static_cast<ostringstream & >( ostringstream() << dec
	//	<< vals->first << EnumDelims[2] << vals->second ).str();
}

int Options::Option::SetTriedFloat(float val, float min, float max)
{
	if (!val && Sign.Is(tOpt::ALLOW0))	min = 0;
	if (val < min || val > max) {
		cerr << optTitle << NameToStr(true) << SepSCl << sValue << setprecision(4) << SPACE << val
			<< " is out of available range [" << min << HPH << max << "]\n";
		return 1;
	}
	NVal = val;
	return 0;
}

int Options::Option::SetEnum(const char* val)
{
	if (!ValRequired())			// boolean template
		NVal = !NVal;			// invers boolean value
	else if (val) {				// non-boolean template, value is defined
		int ind = GetEnumInd(val);
		if (ind < 0)	return PrintWrong(val);
		NVal = ind + MinNVal;	// add possible minimum value as a shift
	}
	// else val==NULL, so use default
	return 0;
}

int Options::Option::SetComb(char* vals)
{
	char* delim;	// a pointer to the first occurrence of delimiter COMMA in vals
	int	ind;

	NVal = 0;	// reset default value
	// run through given values
	for (char* val = vals; true; val = delim + 1) {
		delim = strchr(val, EnumDelims[1]);
		if (delim)	*delim = cNULL;		// temporary cut 'val' to C string with single value
		ind = GetEnumInd(val);
		// set bitwise val, in which running number of each bit (from right side)
		// corresponds to number of detected value
		if (ind >= 0)	NVal = float(int(NVal) ^ (1 << ind));
		if (delim)	*delim = EnumDelims[1];		// restore 'vals' string
		else	break;							// no delimiter: last or single value
	}
	return ind < 0 ? PrintWrong(vals, ind == -2 ? "wrong delimiter in value" : strEmpty) : 0;
}

int Options::Option::SetPair(const char* vals, bool isInt)
{
	if (!vals)	return 0;	// value isn't stated and it's optional (that was checked before)
	const char* delim = strchr(vals, EnumDelims[2]);	// a pointer to the delimiter ';' in vals

	if (!delim)	return PrintWrong(vals, "missed '" + string(1, EnumDelims[2]) + "' delimiter in value");
	PairVals& lim = *((PairVals*)SVal);

	if (delim != vals) {		// first value is set
		if (!IsValidFloat(vals, isInt, true)
			|| SetTriedFloat(float(atof(vals)), lim.Values(PairVals::MIN).first, lim.Values(PairVals::MAX).first))
			return 1;
		((pairVal*)SVal)->first = NVal;		// set first PairVals element
	}
	if (*(delim + 1)) {			// second value is set
		if (!IsValidFloat(delim + 1, isInt)
			|| SetTriedFloat(float(atof(delim + 1)), lim.Values(PairVals::MIN).second, lim.Values(PairVals::MAX).second))
			return 1;
		((pairVal*)SVal)->second = NVal;		// set first PairVals element
	}
	return 0;
}

string Options::Option::ToStr(bool prVal) const
{
	string res(optTitle);
	res += NameToStr(true);
	if (prVal)	res += sSPACE + string(SVal);
	return res;
}

void Options::Option::Print(bool descr) const
{
	if (Sign.Is(tOpt::HIDDEN))	return;

	const BYTE	TabLEN = 8;
	const bool	fixValType = ValType == tENUM || ValType == tCOMB;
	USHORT	len = 0;		// first len is used as counter of printed chars

	if (descr)	cout << SPACE, len++;	// full way: double blank first
	// *** signature
	{
		const string name = NameToStr(false);
		cout << SPACE << name;
		len += short(1 + name.length());
	}
	// *** option value type
	cout << SPACE, len++;
	if (IsValEsc())	cout << '[', len++;
	if (fixValType)
		len += PrintEnumVals();		// print enum values
	else {
		const char* buffer = TypeNames[ValType];
		if (buffer)
			cout << buffer,							// print value type
			len += USHORT(strlen(buffer));
	}
	if (IsValEsc())	cout << ']', len++;

	// *** description
	if (!descr)	return;
	// align description 
	short cnt = IndentInTabs - len / TabLEN;	// IndentInTabs=3 * TabLEN: right boundary of descriptions
	if (cnt <= 0)	cnt = 1;
	if (len + cnt * TabLEN >= (IndentInTabs + 1) * TabLEN)	// too long option,
		cnt = IndentInTabs, cout << LF;			// description on the next line
	for (; cnt; cnt--)	cout << TAB;

	// print description
	len = USHORT(strlen(Descr));	// from now 'len' is used as length of description string,
									// 'cnt' is used as external enum counter (and cnt = 0)
	{
		string buffer(len + 1, 0);
		PrintSubLine(const_cast<char*>(buffer.c_str()), Descr, strchr(Descr, LF), fixValType ? (const char**)SVal : NULL, &cnt);
	}
	if (AddDescr) {
		if (Descr[len - 1] != LF)	cout << SPACE;
		cout << AddDescr;
	}
	if (Sign.Is(tOpt::OBLIG))	cout << " Required";
	else if (ValType >= tHELP)	cout << " and exit";

	// print default value
	if (ValRequired() && NVal != NO_VAL) {
		switch (ValType) {
		case tENUM:
		case tCOMB:
			if (NVal >= MinNVal)	// do not print default if it is never set by user
				PRINT_IN_PRTHS(((char**)SVal)
					[int(NVal) - int(MinNVal)]);	// offset by min enum value
			break;
		case tPR_INT:
		case tPR_FL:
			PRINT_IN_PRTHS(PairValsToStr((pairVal*)SVal)); break;
		case tCHAR:	PRINT_IN_PRTHS(char(NVal)); break;
		default:	PRINT_IN_PRTHS(NVal);
		}
	}
	else if (SVal != NULL)	PRINT_IN_PRTHS(ValType == tENUM ? "NONE" : SVal);
	cout << LF;
}

BYTE Options::Option::PrintEnumVals() const
{
	if (!ValRequired())	return 0;

	char** vals = (char**)SVal;			// array of val images
	BYTE len = BYTE(strlen(vals[0]));	// number of printed chars
	BYTE cnt = BYTE(MaxNVal);			// number of values
	if (MinNVal)							// is range of values limited from below?
		cnt -= BYTE(MinNVal) - 1;
	cout << '<' << vals[0];			// first val image from array of val images
	for (BYTE i = 1; i < cnt; i++) {
		cout << EnumDelims[ValType - tENUM] << vals[i];
		len += BYTE(strlen(vals[i]));
	}
	cout << '>';
	return len + cnt + 1;	// here cnt denotes the number of printed delimiters
}

int Options::Option::GetEnumInd(const char* val)
{
	int i = 0;

	// check for delimiter
	for (char c = *val; c; c = *(val + ++i))
		if (!isalpha(c))	return -2;
	// detect value
	for (i = 0; i < MaxNVal; i++)
		if (!_stricmp(val, ((const char**)SVal)[i]))
			return i;
	return -1;
}

int Options::Option::PrintWrong(const char* val, const string& msg) const
{
	cerr << ToStr(false) << SepSCl << (msg == strEmpty ? "wrong " + sValue : msg);
	if (val) cerr << SPACE << val;
	cerr << LF;
	return 1;
}

#ifdef DEBUG
void Options::Option::Print() const
{
	cout << Str << TAB;
	if (strlen(Str) < 4)	cout << TAB;
	Sign.Print(); cout << LF;
}
#endif

void Options::Usage::Print(Option* opts) const
{
	if (Opt != NO_DEF)	// output option value
		opts[Opt].Print(false);
	else if (Par) {		// output parameter
		if (IsParOblig)	cout << SPACE << Par;
		else			PRINT_IN_PRTHS(Par);
		if (ParDescr)	// output parameter description
			cout << "\n  " << Par << ": " << ParDescr;
	}
	cout << endl;
}

int Options::CheckObligs()
{
	for (int i = 0; i < OptCount; i++)
		if (List[i].CheckOblig() < 0)	return -1;
	return 1;
}

int Options::SetOption(char* opt, char* val, char* nextItem, int& argInd)
{
	int i, res = 0;
	const bool isWord = opt[0] == HPH;
	char* sopt = opt += isWord;	// sliced option

	// parse possibly united short options
	do {
		for (i = 0; i < OptCount; i++)
			if ((res = List[i].SetVal(sopt, isWord, val, nextItem, argInd)) > 0)
				return res;
			else if (!res)	break;
		if (res < 0)	return PrintAmbigOpt(sopt, isWord, "unknown", opt);
	} while (!isWord && *++sopt);	// next slice of united short options
	return 0;
}

bool Options::Find(const char* opt)
{
	for (int i = 0; i < OptCount; i++)
		if (List[i].Str && !strcmp(List[i].Str, opt))
			return true;
	return false;
}

int Options::PrintAmbigOpt(const char* opt, bool isWord, const char* headMsg, const char* inOpt)
{
	cerr << headMsg << SPACE << optTitle << HPH;
	if (isWord)	cerr << HPH << opt;
	else {
		cerr << *opt;
		if (inOpt && *(inOpt + 1))
			cerr << " in " << HPH << inOpt;	// print only non-signle-char splitted short option
	}
	if (inOpt && Find(inOpt))
		cerr << ". Do you mean " << HPH << HPH << inOpt << '?';
	cerr << LF;
	return 1;
}

int	Options::PrintVersion()
{
	cout << Product::Version
#ifdef _ZLIB
		<< "\tzlib " << ZLIB_VERSION
#endif
		<< endl;
	return 1;
}

int Options::PrintSummary(bool prTitle)
{
	if (prTitle)	cout << Product::Title << ": ";
	cout << Product::Descr << LF;
	return 1;
}

int Options::PrintUsage(bool title)
{
	BYTE i, k;
	if (title)		PrintSummary(true), cout << LF;

	// output 'Usage' section
	cout << "Usage:";
	for (k = 0; k < UsageCount; k++) {
		cout << TAB << Product::Title;	PRINT_IN_PRTHS("options");
		// output available options
		for (i = 0; i < OptCount; i++)
			List[i].PrintOblig();
		// input parameters
		Usages[k].Print(List);
	}
	cout << endl;

	// output options section
	cout << "Options:\n";
	for (k = 0; k < GroupCount; k++) {
		if (OptGroups[k])	cout << OptGroups[k] << ":\n";
		for (i = 0; i < OptCount; i++)
			List[i].PrintGroup(k);
	}
	return int(title);
}

string const Options::CommandLine(int argc, char* argv[])
{
	ostringstream ss;
	int i = 0;
	argc--;
	while (i < argc)	ss << argv[i++] << SPACE;
	ss << argv[i];
	return ss.str();
}

int Options::Parse(int argc, char* argv[], const char* obligPar)
{
	if (argc < 2) { Options::PrintUsage(true); return -1; }		// output tip	
	int i, res = 1;
	char* token, * nextToken;	// option or parameter

	for (i = 1; i < argc; i++) {		// argv[0] is the path to the program
		token = argv[i];
		nextToken = argv[i + 1];
		if (*token != HPH) {			// token is not an option
			if (i < argc - 1 				// not a last token
				&& *nextToken == HPH)		// next token is an option
				cerr << token << ": neither option nor parameter" << LF, res = -1;
			break;
		}
		if (SetOption(
			token + 1,								// option without HYPHEN
			i + 1 <= argc ? nextToken : NULL,		// option's value
			i + 2 > argc ? NULL : argv[i + 2],		// next item or NULL
			i										// current index in argc
		))
		{
			res = -1; break;
		}
	}
	if (res > 0)	res = CheckObligs();			// check required options
	if (res > 0 && obligPar && i == argc) {
		cerr << Missing << obligPar << LF;
		res = -1;
	}
	return i * res;
}

const string Options::GetFileName(int indOpt, const char* defName, const string& ext)
{
	Option opt = List[indOpt];
	if (opt.SVal)
		return string(opt.SVal) + ext;
	if (opt.IsValEsc() && opt.Sign.Is(tOpt::TRIMMED))
		return FS::ShortFileName(FS::FileNameWithoutExt(defName)) + ext;
	return strEmpty;
}

const string Options::GetPartFileName(int opt, const char* defName)
{
	const char* outName = Options::GetSVal(opt);

	if (!outName)
		return FS::FileNameWithoutExt(defName);
	string sOutName = string(outName);
	return FS::IsDirExist(outName) ? FS::MakePath(sOutName) + FS::FileNameWithoutExt(defName) : sOutName;
}

#ifdef DEBUG
void Options::Print()
{
	for (int i = 0; i < OptCount; i++)
	{
		cout << setw(2) << i << "  ";
		List[i].Print();
	}
}
#endif

/************************ end of class Options ************************/

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

const char* FS::CheckedDirName(int opt)
{
	const char* name = Options::GetSVal(opt);
	if (name && !IsDirExist(name))	Err(Err::D_NONE, name).Throw();
	return name;
}

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

void PrintTime(long elapsed, bool watch, bool parentheses, bool isLF)
{
	int hours = elapsed / 60;
	int mins = hours % 60;
	hours /= 60;
	if (parentheses)	dout << '(';
	dout << setfill('0') << right;		// right couse it may be chanched by previuos output
	if (hours)	dout << setw(2) << hours << COLON;			// hours
	if (mins || !watch)	dout << setw(2) << mins << COLON;	// mins
	// secs
	if (watch)	dout << setw(5) << fixed << setprecision(2) << (elapsed - mins * 60);
	else		dout << setw(2) << long(elapsed) % 60;
	if (parentheses)	dout << ')';
	if (isLF)	dout << LF, fflush(stdout);
}

bool	TimerBasic::Enabled = false;

long TimerBasic::GetElapsed() const
{
	time_t stopTime;
	time(&stopTime);
	return (long)difftime(stopTime, _startTime);
}

void TimerBasic::Print(long elapsed, const char* title, bool parentheses, bool isLF)
{
	if (title)	dout << title;
	PrintTime(elapsed, false, parentheses, isLF);
}

/************************  end ofclass TimerBasic ************************/

/************************  class Timer ************************/
clock_t	Timer::_StartCPUClock;

void Timer::Stop(int offset, bool parentheses, bool isLF)
{
	if (_enabled) {
		if (offset)	dout << setw(offset) << SPACE;
		PrintTime(GetElapsed(), NULL, parentheses, isLF);
	}
}

/************************  end of class Timer ************************/
#ifdef _TEST
/************************  class Stopwatch ************************/

void Stopwatch::Stop(const string title) const
{
	if (!_isStarted)		return;
	_sumTime += GetElapsed();
	if (title != strEmpty)	PrintTime(_sumTime, (title + sSPACE).c_str(), false, true);
}

/************************  end of class Stopwatch ************************/
#endif	// _TEST
/************************  class StopwatchCPU ************************/

void StopwatchCPU::Stop(const char* title, bool print, bool isLF)
{
	_sumclock += clock() - _clock;
	if (print) {
		if (title)	dout << title;
		PrintTime(_sumclock / CLOCKS_PER_SEC, false, false, isLF);
	}
}

/************************  end of class StopwatchCPU ************************/

#ifdef _MULTITHREAD
/************************  class Mutex ************************/

bool	Mutex::_active;		// true if the mutex really should work
mutex	Mutex::_mutexes[int(Mutex::eType::NONE)];

/************************  end of class Mutex ************************/
#endif	// _MULTITHREAD

/************************ class Chrom ************************/

const char* Chrom::Abbr = "chr";
const BYTE		Chrom::MaxAbbrNameLength = 5;
#ifndef _FQSTATN
const char* Chrom::Marks = "XYM";
const string	Chrom::UndefName = "UNDEF";
const string	Chrom::Short = "chrom";
const string	Chrom::sTitle = "chromosome";
const BYTE		Chrom::MaxShortNameLength = BYTE(Short.length()) + MaxMarkLength;
const BYTE		Chrom::MaxNamedPosLength = BYTE(strlen(Abbr)) + MaxMarkLength + CHRLEN_CAPAC + 1;
BYTE		Chrom::CustomOpt;

chrid Chrom::cID = UnID;		// user-defined chrom ID
chrid Chrom::firstHeteroID = 0;	// first heterosome (X,Y) ID

chrid Chrom::HeteroID(const char cMark)
{
	if (!IsRelativeID())		return cMark;		// absolute ID
	for (size_t i = 0; i < strlen(Marks); i++)
		if (cMark == Marks[i])
			return chrid(firstHeteroID + i);	// relative ID
	return UnID;
}

chrid Chrom::CaseInsID(const char* cMark)
{
	if (isdigit(*cMark)) {				// autosome
		chrid id = atoi(cMark) - 1;
		return id < firstHeteroID ? id : UnID;
	}
	return CaseInsHeteroID(*cMark);		// heterosome
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
	return isdigit(*(cName += prefixLen)) ? atoi(cName) - 1 : HeteroID(*cName);
}

chrid Chrom::ValidateID(const char* cName, size_t prefixLen)
{
	if (!cName)					return UnID;
	cName += prefixLen;							// skip prefix
	for (int i = 1; i <= MaxMarkLength; i++)
		if (cName[i] == USCORE)	return UnID;	// exclude chroms with '_'

	if (isdigit(*cName)) {						// autosome
		chrid id = atoi(cName);
		if (/*IsRelativeID() && */id > firstHeteroID)	firstHeteroID = id;
		return id - 1;
	}
	return CaseInsHeteroID(*cName);				// heterosome
}

void Chrom::SetCustomID(bool prColon)
{
	const char* mark = Options::GetSVal(CustomOpt);		// null if no chrom is set by user
	if (mark && (cID = CaseInsID(mark)) == UnID) {
		ostringstream ss;
		ss << "no such " << sTitle << " in this genome";
		if (prColon)	ss << SepCl;
		ss << Options::OptionToStr(CustomOpt, true);
		Err(ss.str()).Throw();
	}
}

void Chrom::SetCustomOption(int opt/*, bool absIDNumb*/)
{
	CustomOpt = opt;
	if (!(firstHeteroID/* = !absIDNumb*/))
		cID = ValidateID(Options::GetSVal(opt));	// apply absolute numeration discipline
}

const string AutosomeToStr(chrid cid) { return to_string(cid + 1); }

const string Chrom::Mark(chrid cid)
{
	if (cid == UnID)		return UndefName;
	if (IsRelativeID())
		return cid < firstHeteroID ? AutosomeToStr(cid) :
		(cid > firstHeteroID + 2) ? UndefName : string(1, Marks[cid - firstHeteroID]);
	return cid < '9' ? AutosomeToStr(cid) : to_string(cid);
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

/************************  MemStatus ************************/

//bool	MemStatus::_enable;
//LLONG	MemStatus::_startVolume;
//
//LLONG MemStatus::getAvailMemory()
//{
//#ifdef OS_Windows
//	MEMORYSTATUSEX status;
//	status.dwLength = sizeof(status);
//	GlobalMemoryStatusEx(&status);
//	//return (size_t)status.ullAvailPhys;
//	return (size_t)status.ullAvailVirtual;
//#else
//	long pages = sysconf(_SC_PHYS_PAGES);
//	long page_size = sysconf(_SC_PAGE_SIZE);
//	return pages * page_size;
//#endif
//}
//
//void MemStatus::StartObserve(bool enable)
//{
//	if(_enable = enable)	
//		_startVolume = getAvailMemory(); 
//}
//
//void MemStatus::StopObserve()
//{
//	if(_enable) {
//		size_t rem = getAvailMemory();
//		cout << "Unallocated memory: " << (_startVolume - rem) << LF;
//	}
//}
/************************  end of MemStatus ************************/