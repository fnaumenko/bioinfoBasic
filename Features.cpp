/**********************************************************
Feature.cpp
Last modified: 07/31/2024
***********************************************************/
#include "Features.h"

void Features::Init(
	const char* fName,
	ChromSizes* cSizes,
	BYTE scoreInd,
	bool joinOvrl,
	eOInfo oinfo,
	bool abortInvalid
)
{
	FBedReader file(
		fName,
		cSizes,
		scoreInd,
		joinOvrl ? UniBedReader::eAction::JOIN : UniBedReader::eAction::OMIT,
		oinfo,
		abortInvalid
	);
	size_t estItemCnt = file.EstItemCount();

	if (estItemCnt) {
		ReserveItems(estItemCnt);
		_file = &file;
		file.Pass(*this);
		_file = nullptr;
	}
#ifdef _BIOCC
	_narrowLenDistr = file.NarrowLenDistr();
#endif
	//PrintEst(estItemCnt);
}

void Features::AddChrom(chrid cID, size_t cnt)
{
	if (!cnt)	return;
	const size_t lastInd = _items.size();
	AddVal(cID, ItemIndices(lastInd - cnt, lastInd));
}

bool Features::operator()()
{
	if (_file->IsJoined()) {
		_items.back().End = _file->ItemEnd();
		return false;
	}
#ifdef _FEATR_SCORE
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
	DoForChrItems( cit, [&res](cItemsIter it) { res += it->Length(); } );
	return res;
}

chrlen Features::EnrRegnLength(chrid cID, BYTE multiplier, float fLen) const
{
	cIter it = GetIter(cID);
	return it != cEnd() ? EnrRegnLength(it, multiplier, fLen) : 0;
}

//void Features::ScaleScores()
//{
//	for (auto& c : Container()) {
//		const auto itEnd = ItemsEnd(c.second.Data);
//
//		for (auto it = ItemsBegin(c.second.Data); it != itEnd; it++)
//			it->Value /= _maxScore;		// if score is undef then it become 1
//	}
//}


#ifdef _ISCHIP
chrlen Features::GetMinFeatureLength() const
{
	chrlen len, minLen = CHRLEN_MAX;

	DoForItems([&](cItemsIter it) { if ((len = it->Length()) < minLen) minLen = len; });
	return minLen;
}
#endif

#ifdef _BIOCC
chrlen Features::GetMinDistance() const
{
	chrlen dist, minDist = CHRLEN_MAX;

	for (const auto& c : Container()) {
		const auto itEnd = ItemsEnd(c.second.Data);
		auto it = ItemsBegin(c.second.Data);

		for (chrlen end = it++->End; it != itEnd; end = it++->End)
			if ((dist = it->Start - end) < minDist)	minDist = dist;
	}
	return minDist;
}
#endif

//const chrlen UNDEFINED  = std::numeric_limits<int>::max();
#define UNDEFINED	vUNDEF

bool Features::Expand(chrlen expLen, const ChromSizes* cSizes, UniBedReader::eAction action)
{
	if (!expLen)	return false;
	chrlen	cRmvCnt = 0, tRmvCnt = 0;	// counters of removed items in current chrom and total removed items

	for (auto& c : Container()) {							// loop through chroms
		const chrlen cLen = (cSizes && cSizes->IsFilled()) ? (*cSizes)[c.first] : 0;	// chrom length
		const auto itEnd = ItemsEnd(c.second.Data);
		auto it = ItemsBegin(c.second.Data);

		cRmvCnt = 0;
		for (it++; it != itEnd; it++) {
			it->Expand(expLen, cLen);						// next item: compare to previous
			if (it->Start <= prev(it)->End)					// overlapping feature
				if (action == UniBedReader::eAction::JOIN) {
					cRmvCnt++;
					it->Start = UNDEFINED;					// mark item as removed
					(it - cRmvCnt)->End = it->End;
				}
				else if (action == UniBedReader::eAction::ACCEPT)
					tRmvCnt += cRmvCnt,
					cRmvCnt = 0;
				else if (action == UniBedReader::eAction::ABORT) {
					//Err("overlapping feature with an additional extension of " + to_string(expLen)).Throw(false, true);
					dout << "overlapping feature with an additional expansion of " << expLen << LF;
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
		for (auto& c : Container()) {						// loop through chroms
			ItemIndices& data = c.second.Data;
			const auto itEnd = ItemsEnd(data);
			cRmvCnt = 0;
			for (auto it = ItemsBegin(data); it != itEnd; it++)
				if (it->Start == UNDEFINED)		cRmvCnt++;	// skip removed item
				else			newItems.emplace_back(*it);
			data.FirstInd -= tRmvCnt;						// correct indexes
			data.LastInd -= (tRmvCnt += cRmvCnt);
		}
		_items.swap(newItems);
	}
	return true;
}

void Features::CheckFeaturesLength(chrlen len, const string& lenDef, const char* sender) const
{
	DoForItems([&](cItemsIter it) {
		if (it->Length() < len) {
			ostringstream oss("Feature size ");
			oss << it->Length() << " is less than stated " << lenDef << SPACE << len;
			Err(oss.str(), sender).Throw();
		}
		});
}

//void Features::FillRegions(chrid cID, Regions& regn) const
//{
//	const auto& data = At(cID).Data;
//	const auto itEnd = ItemsEnd(data);
//
//	regn.Reserve(ItemsCount(data));
//	for (auto it = ItemsBegin(data); it != itEnd; it++)
//		regn.Add(*it);
//}
