/**********************************************************
FqReader.cpp
Last modified: 10/31/2024
***********************************************************/

#include "FqReader.h"

readlen FqReader::ReadLength() const
{
	//CheckGettingRecord();
	return readlen(LineLengthByInd(READ));
}

const char* FqReader::GetCurrRead() const
{
	//CheckGettingRecord();
	return RealRecord() + LineLengthByInd(HEADER1);
}

const char* FqReader::GetSequence()
{
	const char* record = GetNextRecord();
	if (record != NULL) {
		if (*record != AT)
			Err("non '@' marker; missed header line", LineNumbToStr()).Throw();
		if (*(record + LineLengthByInd(HEADER1, false) + LineLengthByInd(READ, false)) != PLUS)
			Err("non '+' marker; missed second header line", LineNumbToStr()).Throw();
	}
	return record;
}

#ifdef MY_DEBUG
reclen FqReader::RecordLength() const
{
	return LineLengthByInd(HEADER1) + LineLengthByInd(HEADER2) + 2 * LineLengthByInd(READ)
		+ 3;	// 3*LF
}

void FqReader::Print(UINT recCnt)
{
	const char* seq;
	for (UINT cnt = 0; seq = GetSequence(); cnt++)
		if (!recCnt || cnt < recCnt)
			printf("%.*s\n", RecordLength(), seq);
}
#endif

