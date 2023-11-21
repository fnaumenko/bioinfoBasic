/**********************************************************
Feature.cpp
Last modified: 11/21/2023
***********************************************************/
#include "Feature.h"

void Features::AddChrom(chrid cID, size_t cnt)
{
	if (!cnt)	return;
	const chrlen lastInd = chrlen(_items.size());
	AddVal(cID, ItemIndices(lastInd - chrlen(cnt), lastInd));
}

bool Features::operator()()
{
	if (_file->IsJoined()) {
		_items.back().End = _file->ItemEnd();
		return false;
	}
#ifdef _ISCHIP
	float score;
	if (_uniScore) _maxScore = score = 1;
	else {
		score = _file->ItemValue();
		if (score < 0)	_maxScore = score = _uniScore = 1;		// score is undefined in input data
		else if (score > _maxScore)	_maxScore = score;
	}
	_items.emplace_back(_file->ItemRegion(), score);
#else
	_items.emplace_back(_file->ItemRegion());
#endif
	return true;
}

chrlen Features::FeaturesLength(cIter cit) const
{
	chrlen res = 0;
	ForChrItems(cit, [&res](cItemsIter it) { res += it->Length(); });
	return res;
}

chrlen Features::EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const
{
	cIter it = GetIter(cID);
	return it != cEnd() ? EnrRegnLength(it, multiplier, fLen) : 0;
}

#ifdef _ISCHIP

// Scales defined score through all features to the part of 1.
void Features::ScaleScores()
{
	for (auto& c : Container()) {
		const auto itEnd = ItemsEnd(c.second.Data);

		for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
			it->Value /= _maxScore;		// if score is undef then it become 1
	}
}
#else	// NO _ISCHIP

// Copies feature coordinates to external DefRegions.
void Features::FillRegions(chrid cID, Regions& regn) const
{
	const auto& data = At(cID).Data;
	const auto itEnd = ItemsEnd(data);

	regn.Reserve(ItemsCount(data));
	for (auto it = ItemsBegin(data); it != itEnd; it++)
		regn.Add(*it);
}
#endif	// _ISCHIP

chrlen Features::GetMinFeatureLength() const
{
	chrlen len, minLen = CHRLEN_MAX;

	ForAllItems([&](cItemsIter it) { if ((len = it->Length()) < minLen) minLen = len; });
	return minLen;
}


chrlen Features::GetMinDistance() const
{
	chrlen dist, minDist = CHRLEN_MAX;

	for (const auto& c : Container()) {
		const auto itEnd = ItemsEnd(c.second.Data);	// ref ???
		auto it = ItemsBegin(c.second.Data);

		for (chrlen end = it++->End; it != itEnd; end = it++->End)
			if ((dist = it->Start - end) < minDist)	minDist = dist;
	}
	return minDist;
}

//const chrlen UNDEFINED  = std::numeric_limits<int>::max();
#define UNDEFINED	vUNDEF

bool Features::Extend(chrlen extLen, const ChromSizes& cSizes, UniBedInFile::eAction action)
{
	if (!extLen)	return false;
	chrlen	cRmvCnt = 0, tRmvCnt = 0;	// counters of removed items in current chrom and total removed items

	for (auto& c : Container()) {							// loop through chroms
		const chrlen cLen = cSizes.IsFilled() ? cSizes[c.first] : 0;	// chrom length
		const auto itEnd = ItemsEnd(c.second.Data);
		auto it = ItemsBegin(c.second.Data);

		it->Extend(extLen, cLen);			// first item
		cRmvCnt = 0;
		for (it++; it != itEnd; it++) {
			it->Extend(extLen, cLen);						// next item: compare to previous
			if (it->Start <= prev(it)->End)				// overlapping feature
				if (action == UniBedInFile::eAction::JOIN) {
					cRmvCnt++;
					it->Start = UNDEFINED;					// mark item as removed
					(it - cRmvCnt)->End = it->End;
				}
				else if (action == UniBedInFile::eAction::ACCEPT)
					tRmvCnt += cRmvCnt,
					cRmvCnt = 0;
				else if (action == UniBedInFile::eAction::ABORT) {
					//Err("overlapping feature with an additional extension of " + to_string(extLen)).Throw(false, true);
					dout << "overlapping feature with an additional extension of " << extLen << LF;
					return false;
				}
				else if (prev(it)->Start != UNDEFINED)		// OMIT: unmarked item
					cRmvCnt++,
					it->Start = UNDEFINED;	// mark item as removed
		}
	}
	if (cRmvCnt) {		// get rid of items marked as removed 
		vector<Featr> newItems;
		newItems.reserve(_items.size() - tRmvCnt);
		tRmvCnt = 0;
		for (auto& c : Container()) {							// loop through chroms
			ItemIndices& data = c.second.Data;
			const auto itEnd = ItemsEnd(data);
			cRmvCnt = 0;
			for (auto it = ItemsBegin(data); it != itEnd; it++)
				if (it->Start == UNDEFINED)		cRmvCnt++;	// skip removed item
				else			newItems.emplace_back(*it);
			data.FirstInd -= tRmvCnt;				// correct indexes
			data.LastInd -= (tRmvCnt += cRmvCnt);
		}
		_items.swap(newItems);
	}
	return true;
}

void Features::CheckFeaturesLength(chrlen len, const string& lenDef, const char* sender) const
{
	ForAllItems([&](cItemsIter it) {
		if (it->Length() < len) {
			ostringstream oss("Feature size ");
			oss << it->Length() << " is less than stated " << lenDef << SPACE << len;
			Err(oss.str(), sender).Throw();
		}
		});
}
