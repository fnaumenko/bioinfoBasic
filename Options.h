/**********************************************************
Options.h
Provides managing executable options
Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 05/01/2024
***********************************************************/
#pragma once

#include "common.h"

// 'Options' implements main() options and parameters treatment.
static class Options
{
private:
	// types of option values
	// do not forget to support the correlation with Options::TypeNames []
	// tVERS option can be absent, so this 'special' group should be started by tHELP
	enum valType {	// using 'enum class' in this case is inconvenient
		tUNDEF,		// not used
		tNAME,		// string value
		tCHAR,		// character value
		tINT,		// integer value
		tFLOAT,		// float value
		tLONG,		// long value
		tENUM,		// enum value
		tCOMB,		// combo value
		tPR_INT,	// tPAIR_INT pair of integers value
		tPR_FL,		// tPAIR_FL pair of floats value
		tHELP,		// help value
		tVERS,		// version value
		tSUMM,		// special value used to call a program from another program
	};

	// Option types
	enum class tOpt {
		NONE = 0,		// no flags
		OBLIG = 0x01,	// true if option is obligatory
		FACULT = 0x02,	// true if option's value is optional
		ALLOW0 = 0x04,	// true if option's value allowed be 0 while checking min value
		HIDDEN = 0x08,	// hidden option (not printed in help)
		TRIMMED = 0x10,	// true if option has been processed already
		WORD = 0x20,	// true if option is stated as multi-chars
	};

	// Option signes
	struct Signs {
	private:
		BYTE signs;		// tOpt as BYTE
	public:
		// Initializes by sign
		Signs(tOpt sign) { signs = BYTE(sign); }
		// Returns true if given sign is set
		bool Is(tOpt sign) const { return (signs & BYTE(sign)) != 0; }	// '!= 0' to avoid warning
		// Sets given sign
		void Set(tOpt sign) { signs |= BYTE(sign); }

#ifdef DEBUG
		void Print() const { cout << int(Is(WORD)) << int(Is(TRIMMED)) << int(Is(HIDDEN)) << int(Is(OBLIG)); }
#endif
	};

	// structure 'Option' keeps options attributes.
	struct Option {
		const char	Char;		// option - character
		const char* Str;		// option - string
		Signs		Sign;		// option's signs
		const valType ValType;	// type of value
		const BYTE	OptGroup;	// option's category
		float		NVal;		// default or established numeric|enum value
		const float MinNVal;	// minimal permissible value;
								// for enum should be a first defined value
		const float MaxNVal;	// maximal permissible value or count of enum values
		const char* SVal;		// string value, or pointer to enum values array (should be cast to char**)
								// If enum option hase ValRequired==false, the pointer dousn't may be NULL
		const char* Descr;		// tip string
		const char* AddDescr;	// additional tip string

		// Return true if option value is escapable (not obligatory)
		bool IsValEsc() const { return Sign.Is(tOpt::FACULT); }

		// Sets option value.
		//	@param opt: option
		//	@param isword: true if option is a word, false if option is a char
		//	@param val: value of option
		//	@param nextItem: next token after opt and val, or NULL
		//	@param argInd: the current index in argc; increased by 1 if value is accepted
		//	@returns: 0 if success, -1 if not found, 1 if option or value is wrong
		int	SetVal(const char* opt, bool isword, char* val, char* nextItem, int& argInd);

		// Check option for obligatory.
		//	@returns: -1 if option is obligatory but not stated, otherwise 1
		int CheckOblig() const;

		// Prints option if it's obligatory
		void PrintOblig() const { if (Sign.Is(tOpt::OBLIG)) Print(false); }

		// Prints option if it belongs to a group g
		void PrintGroup(BYTE g) const { if (OptGroup == g) Print(true); }

		// Returns option name and value optionally
		string ToStr(bool prVal) const;

		// Prints option in full or short way.
		//	@param descr: if true, prints in full way: 
		//	signature, description (marks as Required if needed), default value, otherwise signature only
		void Print(bool descr) const;

#ifdef DEBUG
		void Print() const;
#endif

	private:
		static const BYTE IndentInTabs;		// indentation in tabs before printing descriptions of option fields
		static const char EnumDelims[];		// specifies delimiter for enum [0] and combi [1] values

		// Checks digital value representation. Prints 'wrong val' message in case of failure
		//	@param str: defined (no NULL) string  representing digital value
		//	@returns: true if digital value representation is correct
		bool IsValidFloat(const char* str, bool isInt, bool isPair = false);

		// Returns option's signature as a string
		//	@param asPointed: true if returns signature as it was stated by user
		string NameToStr(bool asPointed) const;

		// Returns string represented pair of value's separated by delimiter.
		const string PairValsToStr(const pairVal* vals) const;

		// returns true if option value is required
		bool ValRequired() const { return MinNVal != vUNDEF; }

		// Checks limits and set numerical value
		//	@param val: numerical value
		//	@returns: 1 if limits are exceeded, otherwise 0
		int SetTriedFloat(float val, float min, float max);

		// Checks and sets enum option value. Prints 'wrong val' message in case of failure
		//	@param val: input value as C string or NULL if value is optional and not defined
		//	@returns: 0 if success, 1 if wrong value
		int SetEnum(const char* val);

		// Checks and sets enum option value. Prints 'wrong val' message in case of failure
		//	@param val: input value as C string
		//	@returns: 0 if success, 1 if wrong value
		int SetComb(char* val);

		// Checks and sets pair option value
		//	@param vals: pair of values as C string or NULL if value isn't set
		//	@returns: 0 if success, 1 if wrong values
		int SetPair(const char* vals, bool isInt);

		// Recursively prints string with replaced ENUM_REPLACE symbol by enum/combi value.
		//	@param buff: external buffer to copy and output temporary string
		//	@param vals: enum/combi values or NULL for other types
		//	@param cnt: external counter of enum/combi values
		static void PrintTransformDescr(char* buff, const char** vals, short* cnt);

		// Recursively prints string with LF inside as a set of left-margin strings
		// Used to output aligned option descriptions
		// First string is printed from current stdout position.
		// Last substring doesn't include LF.
		//	@param buff: external buffer to copy and output temporary string
		//	@param str: input string with possible LFs
		//	@param subStr: substring of input string to the first LF, or NULL if input string is not ended by LF
		//	@param vals: enum/combi values or NULL for other types
		//	@param cnt: external counter of enum/combi values
		static void PrintSubLine(char* buff, const char* str, const char* subStr, const char** vals, short* cnt);

		// Prints enum or combi values
		//	@returns: number of printed symbols
		BYTE PrintEnumVals() const;

		// Performs a case-insensitive search of given string value among enum values.
		//	@param val: input value as string
		//	@returns: index of finded value in enum,
		//	or -1 if the value is not present in enum,
		//	or -2 if the wrong delimiter is encountered
		int GetEnumInd(const char* val);

		// Ouptuts option with error message to cerr
		//	@param val: value or NULL
		//	@param msg: error message about value
		//	@returns: always 1
		int PrintWrong(const char* val, const string& msg = strEmpty) const;
	};

	// structure 'Usage' is used to output some Usage variants in PrintUsage()
	struct Usage {
		const int	Opt;		// option that should be printed in Usage as obligatory
		const char* Par;		// prog parameter that should be printed in Usage as obligatory
		const bool	IsParOblig;	// true if parameter is obligatory, otherwise printed in square brackets
		const char* ParDescr;	// description of prog parameter that should be printed in Usage

		// Prints Usage params
		void Print(Option* opts) const;
	};

	static const char* sPrSummary;		// summary string printed in help
	static const char* sPrTime;			// run time string printed in help
	static const char* sPrUsage;		// usage string printed in help
	static const char* sPrVersion;		// version string printed in help
	static const char* TypeNames[];	// names of option value types in help
	static const char* OptGroups[];	// names of option groups in help
	static const Usage	Usages[];	// content of 'Usage' variants in help
	static const BYTE	OptCount,		// count of options
		GroupCount,		// count of option groups in help
		UsageCount;		// count of 'Usage' variants in help
	static 	Option		List[];			// list of options. Option 'help' always shuld be
										// the last one, option 'version' - before last.

	// Check obligatory options and output message about first absent obligatory option.
	//	@returns: -1 if some of obligatory options does not exists, otherwise 1
	static int	CheckObligs();

	// Set option [with value] or splitted short options
	//	@param opt: option without HYPHEN
	//	@param val: option's value
	//	@param nextItem: next token after opt and val, or NULL
	//	@param argInd: the current index in argc; increased by 1 if value is accepted
	//	@returns: 0 - success, 1 - false
	static int	SetOption(char* opt, char* val, char* nextItem, int& argInd);

	// Returns true if long option opt is defined
	static bool Find(const char* opt);

	// Ouptuts ambiguous option with error message to cerr
	//	@param opt: option
	//	@param isWord: true if option is a word
	//	@param headMsg: message at the beginning
	//	@param inOpt: initial option (in case of ambiguous composite)
	//	@returns: always 1
	static int PrintAmbigOpt(const char* opt, bool isWord, const char* headMsg, const char* inOpt = NULL);

	// Prints version
	//	@returns: always 1
	static int PrintVersion();

	static int PrintSummary(bool prTitle);

public:
	// Prints 'usage' information
	//	@param title: if true prints title before information
	//	@returns: 1 if title is settinf to true, 0 otherwise
	static int PrintUsage(bool title);

	// Returns option name [and value]
	static string OptionToStr(int opt, bool prVal = false) { return List[opt].ToStr(prVal); }

	// Returns command line.
	//	@param argc: count of main() parameters
	//	@param argv: array of main() parameters
	static const string CommandLine(int argc, char* argv[]);

	// Reset int option value to 0
	static void ResetIntVal(int opt) { List[opt].NVal = 0; }

	// Parses and checks main() parameters and their values.
	//	DataWriter message if some of them is wrong.
	//	@param argc: count of main() pearmeters
	//	@param argv: array of main() pearmeters
	//	@param obligPar: name of required application parameter or NULL if not required
	//	@returns:  index of first parameter (not option) in argv[],
	//	or argc if it is absent,
	//	or negative if tokenize complets wrong
	static int Parse(int argc, char* argv[], const char* obligPar = NULL);

	// Get float value by index
	static float GetFVal(int opt) { return List[opt].NVal; }
	// Get string value by index
	static const char* GetSVal(int opt) { return List[opt].SVal; }
	// Get booling value by index
	static bool GetBVal(int opt) { return List[opt].NVal != 0; }
	// Get UINT value by index
	static UINT GetUIVal(int opt) { return UINT(List[opt].NVal); }
	// Get int value by index
	static int GetIVal(int opt) { return int(List[opt].NVal); }
#ifdef _READS
	// Get read duplicates level
	static char GetRDuplLevel(int opt) { return GetBVal(opt) ? vUNDEF : 0; }
#endif

	// Returns true if the option value is assigned by user
	static bool Assigned(int opt) { return List[opt].Sign.Is(tOpt::TRIMMED); }
	// True if maximal enum value is setting
	static bool IsMaxEnum(int opt) { return List[opt].NVal == List[opt].MaxNVal - 1; }
	// Get maximal permissible numeric value by index
	static UINT GetMaxIVal(int opt) { return UINT(List[opt].MaxNVal); }
	// Returns C-string 'ON' or 'OFF'
	static const char* BoolToStr(int opt) { return Booleans[GetBVal(opt)]; }

#ifdef DEBUG
	static void Print();
#endif
} options;
