/**********************************************************
Region intersection handling and additional classes
2024 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 07/27/2024
-------------------------
***********************************************************/
#pragma once

#include "common.h"

template<typename T>
void DiscardNonOverlapRegions(T rgns[2], fraglen minOverlapLen)
{
	const typename T::iterator itEnd[2]{ rgns[0].end(), rgns[1].end() };
	typename T::iterator it[2]{ rgns[0].begin(), rgns[1].begin() };
	BYTE s;		// index of the left started region: 0 - direct, 1 - reverse
	BYTE suspended = 0;	// if 1st|2nd bit is set then direct|reverse region is suspended and should be analyzed in the next pass

	auto start = [&it](BYTE s) -> chrlen { return T::Start(it[s]); };
	auto end = [&it](BYTE s) -> chrlen { return T::End(it[s]); };

	// unconditional discarde the region
	auto discardeRgn = [&](BYTE s) {
		rgns[s].Discard(it[s]);
		suspended &= ~(1 << s);	// reset suspended s
		it[s]++;
	};
	// unconditional discarde the regions; always applied to the left region
	auto discardeLastRgns = [&it, &itEnd, &discardeRgn](BYTE s) {
		while (it[s] != itEnd[s])
			discardeRgn(s);
	};

	while (it[0] != itEnd[0] && it[1] != itEnd[1]) {
		s = start(0) > start(1);
		if (end(s) > start(!s) + minOverlapLen) {	// strong intersection
			/***************************************
				-----    ?????	more regions?		left overlaps
			+++++++++++++		suspended...		left overlaps
				---------		suspended...		right overlaps
			+++++++	  ?????		more regions?		right overlaps
			***************************************/
			bool valid = true;
			if (T::IsWeak(it[s]))	discardeRgn(s), valid = false;
			if (T::IsWeak(it[!s]))	discardeRgn(!s), valid = false;
			if (valid) {
				rgns[0].Accept(it);
				s ^= end(s) < end(!s);	// flip by condition
				it[!s]++;
				suspended |= 1 << s;	// set suspended s
			}
		}
		else {							// weak intersection or no one
			// conditional close of the region
			if (suspended)	suspended = 0;
			else			rgns[s].Discard(it[s]);

			// close remaining complementary regions
			if (++it[s] == rgns[s].end()) {
				discardeLastRgns(!s);
				break;			// no need to check in while()
			}
		}
	}

	// close 'out of scope' regions
	if ((s = it[1] != itEnd[1]) || it[0] != itEnd[0]) {
		it[s]++;
		discardeLastRgns(s);
	}
}

// 'IGVlocus' is designed to print locus that can be pasted into the address bar of Integrative Genomics Viewer
class IGVlocus
{
	std::string _chrom;
	mutable char _buf[2 * 10 + 5 + 2]{ 0 };	// 2 * max position length + Chrom::MaxAbbrNameLength + 2 separators

	// Prints IGV locus to inner buffer
	void NPrint(chrlen start, chrlen end) const
	{
		const fraglen ROI_ext = 500;
		sprintf(_buf, "%s:%d-%d", _chrom.c_str(), start - ROI_ext, end + ROI_ext);
	}

public:
	IGVlocus(chrid cID = Chrom::UnID) : _chrom(Chrom::AbbrName(cID)) {}

	void SetChrom(chrid cID) { _chrom = Chrom::AbbrName(cID); }

	const char* ChromAbbrName() const { return _chrom.c_str(); }

	// Prints IGV locus to inner buffer
	//	@param start: feature's start position
	//	@param end: feature's end position
	//	@returns: inner buffer
	const char* Print(chrlen start, chrlen end) const { NPrint(start, end); return _buf; }

	// Prints IGV locus to inner buffer
	//	@param pos: feature's start position
	//	@returns: inner buffer
	const char* Print(chrlen pos) const { return Print(pos, pos); }
};

// 'Formatted Writer' - simple output text file with formatted line writer
class FormWriter
{
	FILE* _file;

public:
	FormWriter(const char* name) { _file = fopen(name, "w"); }

	~FormWriter() { fclose(_file); }

	template<typename... Args>
	void Write(const char* format, Args ... args) { fprintf(_file, format, args...); }
};
