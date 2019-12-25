#pragma once

#include<iostream>
#include<sstream>
#include <unordered_map>
#include "BamWalker.h"
#include "RefGenome.h"
#include <thread>
#include <mutex>
#include<fstream>
#include <algorithm>

#include "ssw_cpp.h"

struct rawBkp
{
	std::string from_ref;
	std::string to_ref;
	std::vector<int> from_bkp_list;
	std::vector<int> to_bkp_list;
	bool reverse;
};

struct accBkp
{
	std::string from_ref;
	std::string to_ref;
	std::string from_bkp;
	std::string to_bkp;
	std::string from_side;
	std::string to_side;
	std::string reverse;
	std::string interval;
	std::string read_str;
	std::string ref_str;
	std::string similarity;
};


class GetAccBkp
{
public:
	GetAccBkp();
	~GetAccBkp();


	//std::shared_ptr<RefGenome> ref_genome;
	static std::mutex g_mutex;

	std::vector<rawBkp> raw_breakpoints_list;
	std::vector<accBkp> acc_breakpoints_list;

	static std::unordered_map <std::string, int> RefSeqLength;


	void loadRefGenomeFai(std::string ref_genome_index_file);

	static std::string getMateSeq(std::string sequence);

	void getRawBreakpoints(std::string raw_bkp_file, std::string acc_bkp_file, int rlen, int insert_size);

	static void getAccurateBreakpoints(int start, int end, std::vector<rawBkp>raw_breakpoints_list, std::unordered_map < std::string, std::unordered_map<int, std::vector < std::string >> > left_clipped_reads_dict,
		std::unordered_map<std::string, std::unordered_map<int, std::vector < std::string >>>right_clipped_reads_dict, std::string output, unsigned int threshold, std::string ref_file, int rlen, int insert_size);
	
	void spawnThreads(int n, std::vector<rawBkp>raw_breakpoints_list, std::unordered_map < std::string, std::unordered_map<int, std::vector < std::string >> > left_clipped_reads_dict,
		std::unordered_map<std::string, std::unordered_map<int, std::vector < std::string >>>right_clipped_reads_dict, std::string output, unsigned int threshold, std::string ref_file, int rlen, int insert_size);

	static int sswAlignment(std::string s1, std::string s2);
};


