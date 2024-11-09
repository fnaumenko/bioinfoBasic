/**********************************************************
OrderedData.cpp
Last modified: 11/09/2024
***********************************************************/

#include "OrderedData.h"
#include <fstream>

/************************ AccumCover ************************/
#ifdef MY_DEBUG
void AccumCover::CheckHint(const covmap::iterator& hint, chrlen newPos) const
{
	if (hint != begin() && prev(hint)->first > newPos)
		printf("check hint: new position %d < prev hint %d\n", newPos, prev(hint)->first);
	else if(hint != end() && hint->first < newPos)
		printf("check hint: new position %d > hint %d\n", newPos, hint->first);
}
#endif

void AccumCover::AddRegion(const Region& rgn, bool incrStart)
{
	if (empty()) {
		emplace_hint(end(), rgn.Start, 1);
		emplace_hint(end(), rgn.End, 0);
		return;
	}

	// *** set up 'end' entry
	covmap::iterator it2;		// 'end' entry iterator
	if (incrStart) {
		for (it2 = prev(end()); rgn.End <= it2->first; it2--);

		it2 = emplace_hint(next(it2), rgn.End, it2->second);	// duplicate doesn't change anything
	}
	else {
		it2 = lower_bound(rgn.End);		// 'end' entry iterator

		if (it2 == begin()) {			// insert the entire fragment at the beginning
			if (it2->first != rgn.End)
				it2 = emplace_hint(it2, rgn.End, 0);
			emplace_hint(it2, rgn.Start, 1);
			return;
		}

		if (it2 == end())
			it2 = emplace_hint(it2, rgn.End, 0);					// new last 'end' entry
		else if (it2->first != rgn.End)
			it2 = emplace_hint(it2, rgn.End, prev(it2)->second);	// new 'end' entry
	}

	// *** set up 'start' entry
	auto it1 = prev(it2);							// 'start' entries iterator
	for (; it1 != end() && rgn.Start < it1->first; it1--);

	if (it1 == end())
		it1 = emplace_hint(begin(), rgn.Start, 1);	// new first 'start' entry
	else if (it1->first != rgn.Start)
		it1 = emplace_hint(next(it1), rgn.Start, it1->second + 1);	// new 'start' entry
	else {
		it1->second++;								// incr val at existed 'start' entry
		auto it = prev(it1);
		if (it->second == it1->second)			// previous and current entries have the same value
			erase(it1), it1 = it;				// remove current entry as duplicated
	}

	// *** correct range between 'start' and 'end', set 'end' entry value
	for (it1++; it1 != it2; it1++)				// correct values within range (except the 'end')
		++it1->second;							// increase value
	if ((--it1)->second == it2->second)			// is the last added entry a duplicate?
		erase(it2);								// remove duplicated entry
}

#ifdef _WIG_READER
void AccumCover::AddNextRegion(const Region& rgn, coval val)
{
	if(empty())
		emplace_hint(end(), rgn.Start, val);
	else {
		auto it = prev(end());
		if (it->first == rgn.Start)
			it->second = val;
		else
			emplace_hint(end(), rgn.Start, val);
	}
	emplace_hint(end(), rgn.End, 0);
}
#endif
#ifdef MY_DEBUG
void AccumCover::Print(chrlen maxPos) const
{
	cout << "pos\tval\n";
	DoWithItem([maxPos](const auto& item) {
		if (maxPos && item.first > maxPos)	return;
		cout << item.first << TAB << item.second << LF;
		}
	);
}

void AccumCover::BgPrint(chrlen maxPos) const
{
	cout << "start\tend\tval\n";
	DoWith2Items([maxPos](const auto& it0, const auto& it1) {
		if (maxPos && it0->first > maxPos)	return;
		cout << it0->first << TAB << it1->first << TAB << it0->second << LF;
		}
	);
}
#endif

/************************ AccumCover: end ************************/

/************************ RegionWriter ************************/

const char* RegionWriter::sGRAY = "Silver";	// "175,175,175";

RegionWriter::RegionWriter(FT::eType ftype, eStrand strand, const TrackFields& fields)
	: TxtWriter(ftype, fields.Name, TAB)
{
	static const char* wigFormats[] { FT::BedGraphTYPE, FT::WigTYPE, FT::WigTYPE };
	static const char* StrandCOLORS[2][3] {
		{ "128,128,128", "197,74,74", "0,118,188" },	// grey, red, blue (50,130,190 - foggy blue)
		{ "51,51,51", "102,0,51", "0,0,102" }			// dark grey, dark red, dark blue
	};
	const reclen bufLen = ftype == FT::BED ?
		1000 :		// to save BS bed with extended feilds
		ftype == FT::WIG_FIX ? 300 : 500;

	SetLineBuff(bufLen);
	if (fields.CommLine)	CommLineToIOBuff(*fields.CommLine);

	ostringstream oss;
	oss << "track";
	if(ftype != FT::BED)
		oss << " type=" << wigFormats[int(ftype) - int(FT::BGRAPH)];

	oss << " name=\"" << FS::ShortFileName(fields.Name) << "\" ";
	if (fields.Descr || strand != TOTAL) {
		oss << "description=\"";
		if (fields.Descr)		oss << fields.Descr;
		if (strand != TOTAL)	oss << SepCl << sStrandTITLES[strand] << " strand";
		oss << "\"";
	}
	const char* color = NULL;

	if (ftype != FT::BED) {
		oss << " autoScale=on";
		color = fields.Color ? fields.Color : StrandCOLORS[fields.Shade][strand];
	}
	else {
		if (fields.UseScore)	oss << " useScore=1";
		if (fields.ItemRgb)		oss << " itemRgb=\"On\"";
		if (fields.Color)		color = fields.Color;
		else if(strand!=TOTAL)	color = StrandCOLORS[LIGHT][strand];
	}
	if(color)	oss << " color=" << color;
	StrToIOBuff(oss.str());
}

reclen RegionWriter::AddChromToLine(chrid cID)
{
	LineSetOffset();
	return LineAddStr(Chrom::AbbrName(cID));
}


/************************ RegionWriter: end ************************/

/************************ WigWriter ************************/

void WigWriter::WriteFixStepDeclLine(chrid cID, chrlen pos)
{
	LineSetOffset();
	StrToIOBuff(FT::WigFixSTEP + ChromMarker(cID) + " start=" + to_string(pos));
}

void WigWriter::WriteChromVarStepData(chrid cID, const covmap& cover)
{
	// write declaration line
	LineSetOffset();
	StrToIOBuff(FT::WigVarSTEP + ChromMarker(cID) + " span=1");

	// write data lines
	for (const auto& f : cover) {
		LineAddInts(f.first, f.second, false);	// pos, frequency
		LineToIOBuff();
	}
}

void WigWriter::WriteFixStepRange(chrid cID, chrlen pos, const vector<float>& vals, bool closure)
{
	WriteFixStepDeclLine(cID, pos - bool(vals.front()));

	if (vals.front())
		LineAddSingleFloat(0);		// add zero value to 'open' the curve for the IGV view
	for (float v : vals)
		LineAddSingleFloat(v);
	if (closure && vals.back())
		LineAddSingleFloat(0);		// add zero value to 'close' the curve for the IGV view
}

void WigWriter::WriteFixStepLine(chrid cID, chrlen pos, chrlen ptCnt, float shift)
{
	const bool reversed = shift > 0;
	float val = reversed ? shift : -shift * ptCnt;

	WriteFixStepDeclLine(cID, pos - !reversed);

	LineAddSingleFloat(0);		// if started with unzero value, add zero value to 'open' the curve for the IGV view
	for (chrlen i = 0; i < ptCnt; val += shift, i++)
		LineAddSingleFloat(val);
	if (reversed)
		LineAddSingleFloat(0);	// add zero value to 'close' the curve for the IGV view
}

/************************ WigWriter: end ************************/

/************************ BedGrWriter ************************/

void BedGrWriter::WriteChromData(chrid cID, const covmap& cover)
{
	if (cover.empty())	return;

	const reclen offset = AddChromToLine(cID);

	auto it0 = cover.cbegin(), it = it0;
	const auto end = cover.cend();

	for (++it; it != end; it0 = it++)
		if (it0->second)
			LineAddUInts(it0->first, it->first, it0->second, false),		// start, end, coverage
			LineToIOBuff(offset);
}

/************************ BedGrWriter: end ************************/
