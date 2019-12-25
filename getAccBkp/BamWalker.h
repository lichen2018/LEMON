#pragma once

#include <string>
#include <map>
#include <cmath>
#include"BamReader.h"

class BamWalker
{
public:
	BamWalker() {};
	BamWalker(std::string unique_bam_file, std::string splitter_bam_file);
	int rlen;
	int insert_size;
	BamReader ubr;
	BamReader sbr;


	std::unordered_map < std::string, std::unordered_map<int, std::vector < std::string >> > left_clipped_reads_dict;
	std::unordered_map<std::string, std::unordered_map<int, std::vector < std::string >>>right_clipped_reads_dict;

	void indexSoftClippedReads();
	void getBamStats();



	~BamWalker();
};

