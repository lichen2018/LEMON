// get_acc_bkp.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "GetAccBkp.h"
#include <getopt.h>
#include "BamReader.h"
#include "BamHeader.h"
#include <string>
#include<thread>
#include "BamWalker.h"
#include <iostream>
#include <fstream>

//char* const l_opt_arg = 'ref_sequence:unique_bam:splitter_bam:raw_bkp:output';
char* const short_options = "r:u:s:b:o:t:";

struct option long_options[] = {
	{ "ref_sequence", 1, NULL, 'r' },
	{ "unique_bam",1, NULL, 'u' },
	{ "splitter_bam", 1, NULL, 's' },
	{ "raw_bkp", 1, NULL, 'b' },
	{ "output", 1, NULL, 'o' },
    { "nthreads", 1, NULL, 't'},
	{ 0, 0, 0, 0 }
};


int main(int argc, char * const argv[])
{
    static std::string ref_file;
    static std::string unique_bam_file;
    static std::string splitter_bam_file;
    static std::string raw_bkp_file;
    static std::string acc_bkp_file;
    static int nthreads;
	//BamReader u;
	//u.Open("/mnt/delta_adam/User/lichen/IBD/H4009/H4009C1.dna.unique.bam");
	//BamReader s;
	char opt;
	while ((opt = getopt_long(argc, argv, short_options, long_options, NULL)) != -1)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (opt)
		{
		case 'r': arg >> ref_file; break;
		case 'u': arg >> unique_bam_file; break;
		case 's': arg >> splitter_bam_file; break;
		case 'b': arg >> raw_bkp_file; break;
		case 'o': arg >> acc_bkp_file; break;
        case 't': arg >> nthreads; break;
		}
	}
	std::ifstream file_exist(acc_bkp_file.c_str());
	if(file_exist)
		return 0;
	unsigned int threshold = 15;
	BamWalker bw = BamWalker(unique_bam_file, splitter_bam_file);

	bw.getBamStats();
	bw.indexSoftClippedReads();
	GetAccBkp worker = GetAccBkp();
	worker.getRawBreakpoints(raw_bkp_file, acc_bkp_file, bw.rlen, bw.insert_size);
	worker.loadRefGenomeFai(ref_file);
	std::cout << bw.rlen << "," << bw.insert_size << "\n";
	std::cout << worker.raw_breakpoints_list.size() << "\n";
	worker.spawnThreads(nthreads, worker.raw_breakpoints_list, bw.left_clipped_reads_dict, bw.right_clipped_reads_dict, acc_bkp_file, 
		threshold, ref_file, bw.rlen, bw.insert_size);

	return 0;

}




