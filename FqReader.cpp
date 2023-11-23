/**********************************************************
FqReader.cpp
Last modified: 11/23/2023
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
	return NextRecord() - RecordLength() + LineLengthByInd(HEADER1);
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

