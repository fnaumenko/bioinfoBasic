/**********************************************************
ChromSeq.cpp  2023 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 11/21/2023
***********************************************************/
#include "ChromSeq.h"

bool ChromSeq::LetGaps = true;	// if true then include gaps at the edges of the ref chrom while reading
bool ChromSeq::StatGaps = false;	// if true sum gaps for statistic output

bool ChromSeq::Init(const string& fName, ChromDefRegions& rgns, bool fill)
{
	_seq = nullptr;
	bool getN = StatGaps || LetGaps || rgns.Empty();	// if true then chrom def regions should be recorded
	FaReader file(fName, rgns.Empty() ? &rgns : nullptr);

	_len = file.ChromLength();
	if (fill) {
		try { _seq = new char[_len]; }
		catch (const bad_alloc&) { Err(Err::F_MEM, fName.c_str()).Throw(); }
		const char* line = file.Line();		// First line is readed by FaReader()
		chrlen lineLen;
		_len = 0;

		do	memcpy(_seq + _len, line, lineLen = file.LineLength()),
			_len += lineLen;
		while (line = file.NextGetLine());
	}
	else if (getN)	while (file.NextGetLine());	// just to fill chrom def regions
	file.CLoseReading();	// only makes sense if chrom def regions were filled
	return getN;
}

ChromSeq::ChromSeq(chrid cID, const ChromSizes& cSizes)
{
	_ID = cID;
	ChromDefRegions rgns(cSizes.ServName(cID));	// read from file or new (empty)

	if (Init(cSizes.RefName(cID) + cSizes.RefExt(), rgns, true) && !rgns.Empty())
		_effDefRgn.Set(rgns.FirstStart(), rgns.LastEnd());
	else
		_effDefRgn.Set(0, Length());
	_gapLen = rgns.GapLen();
}

#if defined _READDENS || defined _BIOCC

ChromSeq::ChromSeq(const string& fName, ChromDefRegions& rgns, short minGapLen)
{
	Init(fName, rgns, false);
	rgns.Combine(minGapLen);
}

#endif

//#if defined _FILE_WRITE && defined DEBUG
//#define FA_LINE_LEN	50	// length of wrtied lines
//
//void ChromSeq::Write(const string & fName, const char *chrName) const
//{
//	FaReader file(fName, chrName);
//	chrlen i, cnt = _len / FA_LINE_LEN;
//	for(i=0; i<cnt; i++)
//		file.AddLine(_seq + i * FA_LINE_LEN, FA_LINE_LEN);
//	file.AddLine(_seq + i * FA_LINE_LEN, _len % FA_LINE_LEN);
//	file.Write();
//}
//#endif	// DEBUG
