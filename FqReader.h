/**********************************************************
FqReader.h
Provides FQ reader functionality
2014 Fedor Naumenko (fedor.naumenko@gmail.com)
Last modified: 11/23/2023
***********************************************************/
#pragma once

#include "TxtFile.h"

// 'FqReader' implements reading file in FQ format.
class FqReader : public TxtReader
{
	enum eLineLenIndex { HEADER1, READ, HEADER2, QUAL };

public:
	// Creates new instance for reading by file name
	FqReader(const string& fileName) : TxtReader(fileName, eAction::READ, 4, false) {}

	// Returns checked length of current readed Read.
	readlen ReadLength() const;

	// Gets checked Read from readed Sequence.
	const char* GetCurrRead() const;

	// Returns checked Sequence.
	const char* GetSequence();

	// Returns count of sequences.
	size_t Count() const { return RecordCount(); }
};
