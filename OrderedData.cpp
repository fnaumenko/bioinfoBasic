/**********************************************************
OrderedData.cpp
Last modified: 05/21/2024
***********************************************************/

#include "OrderedData.h"
#include <fstream>

/************************ AccumCover ************************/

void AccumCover::AddRegion0(const Region& frag)
{
	// version using 'find'
	covmap::iterator it1 = find(frag.Start), it2;	// 'start', 'end' entries iterator

	// *** set up 'start' entry
	if (it1 == end()) {							// 'start' entry doesn't exist
		it2 = it1 = emplace(frag.Start, 1).first;
		if (it1 != begin())						// 'start' entry is not the first one
			it1->second += (--it2)->second;		// correct val by prev point; keep it1 unchanged
	}
	else {
		it1->second++;							// incr val at existed 'start' entry
		if (--(it2 = it1) != end()				// decr it2; previous entry exists
			&& it2->second == it1->second)		// previous and current entries have the same value
			erase(it1), it1 = it2;				// remove current entry as duplicated
	}

	// *** set up 'end' entry
	it2 = find(frag.End);
	fraglen val = 0;							// 'end' entry value
	const bool newEnd = it2 == end();

	if (newEnd) {								// 'end' entry doesn't exist
		it2 = emplace(frag.End, 0).first;
		val = next(it2) == end() ? 1 :			// 'end' entry is the last one
			prev(it2)->second;					// grab val by prev entry
	}

	// *** correct range between 'start' and 'end', set 'end' entry value
	for (it1++; it1 != it2; it1++)				// correct values within range
		val = ++it1->second;					// increase value
	if ((--it1)->second == it2->second)			// is the last added entry a duplicate?
		erase(it2);								// remove duplicated entry
	else if (newEnd)
		it2->second = --val;					// set new 'end' entry value
}
void AccumCover::AddRegion(const Region& frag)
{
	// version using 'lower_bound'
	covmap::iterator it1 = lower_bound(frag.Start), it2;	// 'start', 'end' entries iterator

	// *** set up 'start' entry
	if (it1 == end() || it1->first != frag.Start) {			// 'start' entry doesn't exist
		it2 = it1 = emplace_hint(it1, frag.Start, 1);
		if (it1 != begin())						// 'start' entry is not the first one
			it1->second += (--it2)->second;		// correct val by prev point; keep it1 unchanged
	}
	else {
		it1->second++;							// incr val at existed 'start' entry
		if (--(it2 = it1) != end()				// decr it2; previous entry exists
		&& it2->second == it1->second)			// previous and current entries have the same value
			erase(it1), it1 = it2;				// remove current entry as duplicated
	}

	// *** set up 'end' entry
	it2 = lower_bound(frag.End);
	const bool newEnd = it2 == end() || it2->first != frag.End;
	fraglen val = 0;							// 'end' entry value

	if (newEnd) {								// 'end' entry doesn't exist
		it2 = emplace_hint(it2, frag.End, 0);
		val = next(it2) == end() ? 1 :			// 'end' entry is the last one
			prev(it2)->second;					// grab val by prev entry
	}

	// *** correct range between 'start' and 'end', set 'end' entry value
	for (it1++; it1 != it2; it1++)				// correct values within range
		val = ++it1->second;					// increase value
	if ((--it1)->second == it2->second)			// is the last added entry a duplicate?
		erase(it2);								// remove duplicated entry
	else if (newEnd)
		it2->second = --val;					// set new 'end' entry value
}

#ifdef _WIG
void AccumCover::AddNextRegion(const Region& rgn, coval val)
{
	//if(empty())
	//	emplace_hint(end(), rgn.Start, val);
	//else {
	//	auto it = prev(end());
	//	if (it->first == rgn.Start)
	//		it->second = val;
	//	else
	//		emplace_hint(end(), rgn.Start, val);
	//}
	if (!empty() && prev(end())->first == rgn.Start)
		prev(end())->second = val;
	else
		emplace_hint(end(), rgn.Start, val);
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

//const char* RegionWriter::sTRACK = "track";
const char* RegionWriter::sGRAY = "Silver";	// "175,175,175";

RegionWriter::RegionWriter(FT::eType ftype, eStrand strand, const TrackFields& fields)
	: TxtWriter(ftype, fields.Name, TAB)
{
	static const char* wigFormats[] = { FT::BedGraphTYPE, FT::WigTYPE, FT::WigTYPE };
	static const char* StrandCOLORS[] = { "128,128,128", "197,74,74", "0,118,188" };	// grey, red, blue; 50,130,190 - foggy blue
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
		color = fields.Color ? fields.Color : StrandCOLORS[strand];
	}
	else {
		if (fields.UseScore)	oss << " useScore=1";
		if (fields.ItemRgb)		oss << " itemRgb=\"On\"";
		if (fields.Color)		color = fields.Color;
		else if(strand!=TOTAL)	color = StrandCOLORS[strand];
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
