#include "BamWalker.h"


BamWalker::BamWalker(std::string unique_bam_file, std::string splitter_bam_file)
{
	ubr.Open(unique_bam_file);
	sbr.Open(splitter_bam_file);
}


BamWalker::~BamWalker()
{
}


void BamWalker::getBamStats()
{
	BamRecord r;
	int count = 0;
	int insert_size_sum = 0;
	int read_length_sum = 0;
	std::vector<int>insert_size_list;
	while (ubr.GetNextRecord(r))
	{
		if (r.ProperPair() && r.PairedFlag() && r.InsertSize() > 0 && r.InsertSize() < 1000 && !r.DuplicateFlag() && r.MappedFlag() && r.MateMappedFlag())
		{
			count++;
			read_length_sum += r.Sequence().size();
			insert_size_list.push_back(r.InsertSize());
			insert_size_sum += r.InsertSize();
		}
	}
	rlen = read_length_sum / count;
	double mean = (double)insert_size_sum / count;
	double sdev = 0;
	for (int i = 0; static_cast<unsigned int>(i) < insert_size_list.size(); i++)
	{
		sdev += (insert_size_list[i] - mean)*(insert_size_list[i] - mean);
	}
	sdev = (double)sdev / (count - 1);
	sdev = std::sqrt(sdev);
	insert_size = int(mean + sdev);
}




void BamWalker::indexSoftClippedReads()
{
	BamRecord r;
	while (sbr.GetNextRecord(r))
	{
		std::string reference_name = r.ChrName(sbr.Header());

		if (static_cast<unsigned int>(r.AlignmentEndPosition()) < r.Sequence().size())
		{//right clipped reads
			int bkp_pos = r.Position() + r.AlignmentEndPosition() - r.AlignmentPosition();
			std::string clipped_seq = r.Sequence().substr(r.AlignmentEndPosition(), r.Sequence().size());
			if (right_clipped_reads_dict.find(reference_name) == right_clipped_reads_dict.end())
			{//not find reference name
				std::unordered_map<int, std::vector < std::string >>read_pos_dict;
				std::vector<std::string>reads_vec;
				reads_vec.push_back(clipped_seq);
				read_pos_dict.insert(std::pair<int, std::vector<std::string>>(bkp_pos, reads_vec));
				right_clipped_reads_dict.insert(std::pair < std::string, std::unordered_map<int, std::vector < std::string >>>(reference_name, read_pos_dict));
			}
			else
			{
				if (right_clipped_reads_dict[reference_name].find(bkp_pos) == right_clipped_reads_dict[reference_name].end())
				{
					std::vector<std::string>reads_vec;
					reads_vec.push_back(clipped_seq);
					right_clipped_reads_dict[reference_name].insert(std::pair<int, std::vector<std::string>>(bkp_pos, reads_vec));
				}
				else
				{
					right_clipped_reads_dict[reference_name][bkp_pos].push_back(clipped_seq);
				}
			}
		}
		if (r.AlignmentPosition() > 0)
		{
			int bkp_pos = r.Position();
			std::string clipped_seq = r.Sequence().substr(0, r.AlignmentPosition());
			if (left_clipped_reads_dict.find(reference_name) == left_clipped_reads_dict.end())
			{//not find reference name
				std::unordered_map<int, std::vector < std::string >>read_pos_dict;
				std::vector<std::string>reads_vec;
				reads_vec.push_back(clipped_seq);
				read_pos_dict.insert(std::pair<int, std::vector<std::string>>(bkp_pos, reads_vec));
				left_clipped_reads_dict.insert(std::pair < std::string, std::unordered_map<int, std::vector < std::string >>>(reference_name, read_pos_dict));
			}
			else
			{
				if (left_clipped_reads_dict[reference_name].find(bkp_pos) == left_clipped_reads_dict[reference_name].end())
				{
					std::vector<std::string>reads_vec;
					reads_vec.push_back(clipped_seq);
					left_clipped_reads_dict[reference_name].insert(std::pair<int, std::vector<std::string>>(bkp_pos, reads_vec));
				}
				else
				{
					left_clipped_reads_dict[reference_name][bkp_pos].push_back(clipped_seq);
				}
			}
		}
	}
}