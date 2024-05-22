/**********************************************************
ChromSeq.h  2023 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 05/22/2024
***********************************************************/
#pragma once

#include "ChromData.h"

// 'ChromSeq' represented chromosome as an array of nucleotides
class ChromSeq
{
private:
	chrid	_ID;			// ID of chromosome
	chrlen	_len;			// length of chromosome
	chrlen	_gapLen = 0;	// total length of gaps
	Region	_effDefRgn;		// effective defined region (except 'N' at the begining and at the end)
	char* _seq = nullptr;	// the nucleotides buffer

	// Initializes instance and/or chrom's defined regions
	//	@param fName: file name
	//	@param rgns: chrom's defined regions: ripe or new
	//	@param fill: if true fill sequence and def regions, otherwise def regions only
	//	@returns: true if chrom def regions are stated
	bool Init(const string& fName, ChromDefRegions& rgns, bool fill);

public:
	static bool	LetGaps;	// if true then include gaps at the edges of the ref chrom while reading
	static bool	StatGaps;	// if true count sum gaps for statistic output

	~ChromSeq() { if(_seq)	delete[] _seq; }

	// Gets chrom legth
	chrlen Length()	const { return _len; }

	// Gets Read on position or NULL if exceeding the chrom length
	//const char* Read(chrlen pos) const { return pos > _len ? NULL : _seq + pos;	}

	//const char* Read(chrlen pos, readlen len) const { return pos + len > Length() ? NULL : _seq + pos; }

	// Gets subsequence without exceeding checking 
	const char* Seq(chrlen pos) const { return _seq + pos; }

	// Creates a stub instance (for sampling cutting)
	//	@param len: chrom length
	ChromSeq(chrlen len) : _ID(Chrom::UnID), _seq(NULL), _len(len) { _effDefRgn.Set(0, len); }

	// Creates and fills new instance
	ChromSeq(chrid cID, const ChromSizes& cSizes);

	// Returns chrom ID
	chrid ID() const { return _ID; }

	// Gets count of nucleotides outside of defined region
	chrlen UndefLength() const { return Length() - _effDefRgn.Length(); }

	//float UndefLengthInPerc() const { return 100.f * (Length() - _effDefRgn.Length()) / Length(); }

	// Gets defined region (without N at the begining and at the end)
	const Region& DefRegion() const { return _effDefRgn; }

	// Gets start position of first defined nucleotide
	chrlen Start()	const { return _effDefRgn.Start; }

	// Gets end position of last defined nucleotide
	chrlen End()	const { return _effDefRgn.End; }

	// Gets total length of gaps
	chrlen GapLen()	const { return _gapLen; }

#if defined _READDENS || defined _BIOCC

	// Creates an empty instance and fills chrom's defined regions
	//	@param fName: FA file name with extension
	//	@param rgns: new chrom's defined regions
	//	@param minGapLen: minimal length which defines gap as a real gap
	ChromSeq(const string& fName, ChromDefRegions& rgns, short minGapLen);

#endif

//#if defined _TXT_WRITER && defined DEBUG 
//	// Saves instance to file by fname
//	void Write(const string & fname, const char *chrName) const;
//#endif
};