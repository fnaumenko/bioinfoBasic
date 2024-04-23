/**********************************************************
ChromData.cpp
Last modified: 04/20/2024
***********************************************************/

#include "ChromData.h"
#include <algorithm>    // std::sort
#include <fstream>		// to write ChromSizes without defined _FILE_WRITE


inline int ChromSizes::CommonPrefixLength(const string& fName, BYTE extLen)
{
	// a short file name without extention
	return Chrom::PrefixLength(fName.substr(0, fName.length() - extLen).c_str());
}

void ChromSizes::Read(const string& fName)
{
	TabReader file(fName, FT::CSIZE);	// file check already done

	while (file.GetNextLine()) {
		chrid cID = Chrom::ValidateIDbyAbbrName(file.StrField(0));
		if (cID != Chrom::UnID)
			AddValue(cID, ChromSize(file.LongField(1)));
	}
}

void ChromSizes::Write(const string& fName) const
{
	ofstream file;

	file.open(fName.c_str(), ios_base::out);
	for (cIter it = cBegin(); it != cEnd(); it++)
		file << Chrom::AbbrName(CID(it)) << TAB << Length(it) << LF;
	file.close();
}

chrid ChromSizes::GetChromIDs(vector<chrid>& cIDs, const string& gName)
{
	vector<string> files;
	if (!FS::GetFiles(files, gName, _ext))		return 0;

	chrid	cid;				// chrom ID relevant to current file in files
	int		prefixLen;			// length of prefix of chrom file name
	chrid	extLen = BYTE(_ext.length());
	chrid	cnt = chrid(files.size());

	cIDs.reserve(cnt);
	sort(files.begin(), files.end());
	// remove additional names and sort listFiles
	for (chrid i = 0; i < cnt; i++) {
		if ((prefixLen = CommonPrefixLength(files[i], extLen)) < 0)		// right chrom file name
			continue;
		// filter additional names
		cid = Chrom::ValidateID(files[i].substr(prefixLen, files[i].length() - prefixLen - extLen));
		if (cid != Chrom::UnID) 		// "pure" chrom's name
			cIDs.push_back(cid);
	}
	sort(cIDs.begin(), cIDs.end());
	return chrid(cIDs.size());
}

void ChromSizes::SetPath(const string& gPath, const char* sPath, bool prMsg)
{
	_gPath = FS::MakePath(gPath);
	if (sPath && !FS::CheckDirExist(sPath, false))
		_sPath = FS::MakePath(sPath);
	else
		if (FS::IsDirWritable(_gPath.c_str()))
			_sPath = _gPath;
		else {
			_sPath = strEmpty;
			if (prMsg)
				Err("reference folder closed for writing and service folder is not pointed.\n").
				Warning("Service files will not be saved!");
		}
}

void  ChromSizes::SetTreatedChrom(chrid cID)
{
	if (cID != Chrom::UnID)
		for (auto& c : Chroms::Container())
			c.second.Treated = c.first == cID;
}

ChromSizes::ChromSizes(const char* gName, bool prMsg, const char* sPath, bool checkGRef)
{
	_ext = _gPath = _sPath = strEmpty;

	Chrom::SetRelativeMode();
	if (gName) {
		if (FS::IsDirExist(FS::CheckedFileDirName(gName))) {	// gName is a directory
			_ext = FT::Ext(FT::FA);
			SetPath(gName, sPath, prMsg);
			const string cName = _sPath + FS::LastDirName(gName) + FT::Ext(FT::CSIZE);
			const bool isExist = FS::IsFileExist(cName.c_str());
			vector<chrid> cIDs;		// chrom's ID fill list

			// fill list with inizialised chrom ID and set _ext
			if (!GetChromIDs(cIDs, gName)) {				// fill list from *.fa
				_ext += ZipFileExt;			// if chrom.sizes exists, get out - we don't need a list
				if (!isExist && !GetChromIDs(cIDs, gName))	// fill list from *.fa.gz
					Err(Err::MsgNoFile("*", true, FT::Ext(FT::FA)), gName).Throw();
			}

			if (isExist)	Read(cName);
			else {							// generate chrom.sizes
				for (chrid cid : cIDs)
					AddValue(cid, ChromSize(FaReader(RefName(cid) + _ext).ChromLength()));
				if (IsServAvail())	Write(cName);
				if (prMsg)
					dout << FS::ShortFileName(cName) << SPACE << (IsServAvail() ? "created" : "generated") << LF,
					fflush(stdout);			// std::endl is unacceptable
			}
		}
		else {
			if (checkGRef)	Err("is not a directory", gName).Throw();
			Read(gName);		// gName is a chrom.sizes file
			_sPath = FS::DirName(gName, true);
		}
		Chrom::SetUserCID();

		SetTreatedChrom(Chrom::UserCID());
	}
	else if (sPath)
		_gPath = _sPath = FS::MakePath(sPath);	// initialized be service dir; _ext is empty!
	// else instance remains empty
}

void ChromSizes::Init(const string& headerSAM)
{
	if (!IsFilled())
		Chrom::ValidateIDs(
			headerSAM,
			[this](chrid cID, const char* header) { AddValue(cID, atol(header)); },
			true
		);
}

genlen ChromSizes::GenSize() const
{
	if (!_gsize)
		for (const auto& c : *this)	_gsize += c.second.Data.Real;
	return _gsize;
}

#ifdef MY_DEBUG
void ChromSizes::Print() const
{
	//	cout << "ChromSizes: count: " << int(ChromCount()) << endl;
	//	cout << "ID\tchrom\t";
	////#ifdef _ISCHIP
	////		cout << "autosome\t";
	////#endif
	//		cout << "size\n";
	//	for(cIter it=cBegin(); it!=cEnd(); it++) {
	//		cout << int(CID(it)) << TAB << Chrom::AbbrName(CID(it)) << TAB;
	////#ifdef _ISCHIP
	////		cout << int(IsAutosome(CID(it))) << TAB;
	////#endif
	//		cout << Length(it) << TAB << LF;
	//	}

	printf("ChromSizes: count: %d\n", ChromCount());
	printf("ID  chrom   size      treated\n");
	for (const auto& c : *this)
		printf("%2d  %-8s%9d  %d\n", c.first, Chrom::AbbrName(c.first).c_str(), c.second.Data.Real, c.second.Treated);
}
#endif	// MY_DEBUG
