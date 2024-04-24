#include "Options.h"

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
	if (opt.IsValEsc() /*&& opt.Sign.Is(tOpt::TRIMMED)*/)
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
