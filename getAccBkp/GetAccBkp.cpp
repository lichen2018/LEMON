#include "GetAccBkp.h"
#include <time.h>
#include <thread>
#include <chrono>


std::mutex GetAccBkp:: g_mutex;
std::unordered_map <std::string, int> GetAccBkp::RefSeqLength;
double time_limit = 1500;

GetAccBkp::GetAccBkp()
{
}


GetAccBkp::~GetAccBkp()
{
}



void GetAccBkp::loadRefGenomeFai(std::string ref_genome_index_file)
{
	std::string fai_file = ref_genome_index_file + ".fai";
	std::ifstream file(fai_file);
	std::string line;
	while (std::getline(file, line))
	{
		std::istringstream iss(line);
		std::string str;
		std::string ref_name;
		int length = 0;
		int index = 0;
		while (std::getline(iss, str, '\t'))
		{
			std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
			str.erase(end_pos, str.end());
			if (index == 0)
				ref_name = str;
			if (index == 1)
			{
				length = std::stoi(str);
				if (!(RefSeqLength.find(ref_name) == RefSeqLength.end()))
				{
					RefSeqLength[ref_name] = length;
				}
				break;
			}
			index++;
		}
	}
}


std::string GetAccBkp::getMateSeq(std::string sequence)
{
	std::string mate_seq = "";
	for (int i = sequence.size() - 1; i >= 0; i--)
	{
		if (sequence[i] == 'A')
			mate_seq += 'T';
		else if (sequence[i] == 'T')
			mate_seq += 'A';
		else if (sequence[i] == 'G')
			mate_seq += 'C';
		else if(sequence[i] == 'C')
			mate_seq += 'G';
	}
	return mate_seq;
}

void GetAccBkp::getRawBreakpoints(std::string raw_bkp_file, std::string acc_bkp_file, int rlen, int insert_size)
{
	std::ifstream file;
	file.open(raw_bkp_file, std::ifstream::in);
	std::string line;
	while (std::getline(file, line))
	{
		std::istringstream iss(line);
		std::string str;
		std::string from_ref;
		std::string to_ref;
		int from_pos = 0;
		int from_left_pos = 0;
		int from_right_pos = 0;
		int to_pos = 0;
		int to_left_pos = 0;
		int to_right_pos = 0;
		int index = 0;
		bool reverse = false;
		std::string from_side;
		std::string to_side;
		std::string reverse_flag;
		std::string s1;
		std::string s2;
		std::string tmp_reverse("True");
		std::string read_name;
		std::string interval;
		while (std::getline(iss, str, ','))
		{
			std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
			str.erase(end_pos, str.end());
			if (index == 0)
				from_ref = str;
			if (index == 1)
				from_pos = std::stoi(str);
			if (index == 2)
				from_left_pos = std::stoi(str);
			if (index == 3)
				from_right_pos = std::stoi(str);
			if (index == 4)
				to_ref = str;
			if (index == 5)
				to_pos = std::stoi(str);
			if (index == 6)
				to_left_pos = std::stoi(str);
			if (index == 7)
				to_right_pos = std::stoi(str);
			if (index == 9)
			{
				reverse_flag = str;
				if (str.compare(tmp_reverse) == 0)
					reverse = true;
			}

			index++;
		}
		if (RefSeqLength.find(from_ref) == RefSeqLength.end())
		{
			RefSeqLength.insert(std::pair<std::string, int>(from_ref, 0));
		}
		if (RefSeqLength.find(to_ref) == RefSeqLength.end())
		{
			RefSeqLength.insert(std::pair<std::string, int>(to_ref, 0));
		}

		bool found = false;
		bool f1 = false;
		bool f2 = false;
		index = -1;
		for (int i = 0; static_cast<unsigned int>(i) < raw_breakpoints_list.size(); i++)
		{
			rawBkp tmp = raw_breakpoints_list[i];
			if (tmp.from_ref == from_ref && tmp.to_ref == to_ref && tmp.reverse == reverse)
			{
				for (int j = 0; static_cast<unsigned int>(j) < tmp.from_bkp_list.size(); j++)
				{
					if (f1)
						break;
					for (int k = 0; static_cast<unsigned int>(k) < tmp.to_bkp_list.size(); k++)
					{
						if (std::abs(tmp.from_bkp_list[j] - from_pos) < rlen && std::abs(tmp.to_bkp_list[k] - to_pos) < rlen)
						{
							f1 = true;
							index = i;
							break;
						}
					}
				}
			}
			if (tmp.from_ref == to_ref && tmp.to_ref == from_ref && tmp.reverse == reverse)
			{
				for (int k = 0; static_cast<unsigned int>(k) < tmp.to_bkp_list.size(); k++)
				{
					if (f2)
						break;
					for (int j = 0; static_cast<unsigned int>(j) < tmp.from_bkp_list.size(); j++)
					{
						if (std::abs(tmp.from_bkp_list[j] - to_pos) < rlen && std::abs(tmp.to_bkp_list[k] - from_pos) < rlen)
						{
							f2 = true;
							index = i;
							break;
						}
					}
				}
			}
			if (f1 && f2)
				break;
		}
		if (f1)
		{
			bool found = false;
			for (int i = 0; static_cast<unsigned int>(i) < raw_breakpoints_list[index].from_bkp_list.size(); i++)
			{
				if (raw_breakpoints_list[index].from_bkp_list[i] == from_pos)
				{
					found = true;
					break;
				}
			}
			if (!found)
				raw_breakpoints_list[index].from_bkp_list.push_back(from_pos);
			found = false;
			for (int i = 0; static_cast<unsigned int>(i) < raw_breakpoints_list[index].to_bkp_list.size(); i++)
			{
				if (raw_breakpoints_list[index].to_bkp_list[i] == to_pos)
				{
					found = true;
					break;
				}
			}
			if (!found)
				raw_breakpoints_list[index].to_bkp_list.push_back(to_pos);
			continue;
		}
		if (f2)
		{
			bool found = false;
			for (int i = 0; static_cast<unsigned int>(i) < raw_breakpoints_list[index].from_bkp_list.size(); i++)
			{
				if (raw_breakpoints_list[index].from_bkp_list[i] == to_pos)
				{
					found = true;
					break;
				}
			}
			if (!found)
				raw_breakpoints_list[index].from_bkp_list.push_back(to_pos);
			found = false;
			for (int i = 0; static_cast<unsigned int>(i) < raw_breakpoints_list[index].to_bkp_list.size(); i++)
			{
				if (raw_breakpoints_list[index].to_bkp_list[i] == from_pos)
				{
					found = true;
					break;
				}
			}
			if (!found)
				raw_breakpoints_list[index].to_bkp_list.push_back(from_pos);
			continue;
		}
		if (!f1 && !f2)
		{
			rawBkp tmp_bkp;
			tmp_bkp.from_ref = from_ref;
			tmp_bkp.to_ref = to_ref;
			tmp_bkp.from_bkp_list.push_back(from_pos);
			tmp_bkp.to_bkp_list.push_back(to_pos);
			tmp_bkp.reverse = reverse;
			raw_breakpoints_list.push_back(tmp_bkp);
		}
	}
}


int GetAccBkp::sswAlignment(std::string s1, std::string s2)
{
	int32_t maskLen = strlen(s2.c_str()) / 2;
	maskLen = maskLen < 15 ? 15 : maskLen;
	StripedSmithWaterman::Aligner aligner;

	// Declares a default filter

	StripedSmithWaterman::Filter filter;

	// Declares an alignment that stores the result

	StripedSmithWaterman::Alignment alignment;

	// Aligns the query to the ref
	aligner.Align(s1.c_str(), s2.c_str(), s2.size(), filter, &alignment, maskLen);
	int matches = alignment.matches;

	return matches;
}


void GetAccBkp::getAccurateBreakpoints(int start, int end, std::vector<rawBkp>raw_breakpoints_list, std::unordered_map < std::string, std::unordered_map<int, std::vector < std::string >> > left_clipped_reads_dict,
	std::unordered_map<std::string, std::unordered_map<int, std::vector < std::string >>>right_clipped_reads_dict, std::string output, unsigned int threshold, std::string ref_file, int rlen, int insert_size)
{
	RefGenome* rg = new RefGenome();
	rg->LoadIndex(ref_file);
	std::vector<accBkp>better_breakpoint_list;
	for (int index = start; index < end; index++)
	{
		unsigned int max_match_split_length = 0;
		bool found_flag = false;
		bool reverse = raw_breakpoints_list[index].reverse;
		std::string from_side;
		std::string to_side;
		std::string read_str;
		std::string ref_str;
		double similarity;
		int from_pos;
		int to_pos;
		auto startTime = std::chrono::system_clock::now();
		for (int count = 0; count < 2; count++)
		{
			std::string ref_1;
			std::string ref_2;
			std::vector<int> pos_list_1;
			std::vector<int> pos_list_2;
			auto curTime = std::chrono::system_clock::now();
			std::chrono::duration<double> elapsed_seconds = curTime - startTime;
			if (elapsed_seconds.count() > time_limit && !found_flag)
			{
				break;
			}
			if (found_flag)
				break;
			if (count == 0)
			{
				ref_1 = raw_breakpoints_list[index].from_ref;
				ref_2 = raw_breakpoints_list[index].to_ref;
				pos_list_1 = raw_breakpoints_list[index].from_bkp_list;
				pos_list_2 = raw_breakpoints_list[index].to_bkp_list;
			}
			else
			{
				ref_1 = raw_breakpoints_list[index].to_ref;
				ref_2 = raw_breakpoints_list[index].from_ref;
				pos_list_1 = raw_breakpoints_list[index].to_bkp_list;
				pos_list_2 = raw_breakpoints_list[index].from_bkp_list;
			}
			if (ref_1 == ref_2)
				continue;
			int from_left_low = std::max(0, *std::min_element(pos_list_1.begin(), pos_list_1.end()) - insert_size);
			int from_right_up = *std::max_element(pos_list_1.begin(), pos_list_1.end()) + rlen + insert_size;
			int from_left_up = *std::max_element(pos_list_1.begin(), pos_list_1.end()) + threshold;
			int from_right_low = *std::min_element(pos_list_1.begin(), pos_list_1.end()) - threshold;
			std::vector<int>from_low_bound;
			std::vector<int>from_up_bound;
			from_low_bound.push_back(from_left_low);
			from_low_bound.push_back(from_right_low);
			from_up_bound.push_back(from_left_up);
			from_up_bound.push_back(from_right_up);
			int to_left_low = std::max(0, *std::min_element(pos_list_2.begin(), pos_list_2.end()) - insert_size);
			int to_right_up = *std::max_element(pos_list_2.begin(), pos_list_2.end()) + rlen + insert_size;
			int to_left_up = *std::max_element(pos_list_2.begin(), pos_list_2.end()) + threshold;
			int to_right_low = *std::min_element(pos_list_2.begin(), pos_list_2.end()) - threshold;
			std::vector<int>to_low_bound;
			std::vector<int>to_up_bound;
			to_low_bound.push_back(to_left_low);
			to_low_bound.push_back(to_right_low);
			to_up_bound.push_back(to_left_up);
			to_up_bound.push_back(to_right_up);
			int from_index;
			for (int to_index = 0; to_index < 2; to_index++)
			{
				auto curTime = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = curTime - startTime;
				if (elapsed_seconds.count() > time_limit && !found_flag)
				{
					break;
				}

				if (found_flag)
					break;
				if (!reverse)
				{
					if (to_index == 0)
						from_index = 1;
					else
						from_index = 0;
				}
				else
				{
					if (to_index == 0)
						from_index = 0;
					else
						from_index = 1;
				}
				if (from_index == 0)
				{
					auto curTime = std::chrono::system_clock::now();
					std::chrono::duration<double> elapsed_seconds = curTime - startTime;
					if (elapsed_seconds.count() > time_limit && !found_flag)
					{
						break;
					}
					if (found_flag)
						break;
					if (!(left_clipped_reads_dict.find(ref_1) == left_clipped_reads_dict.end()))
					{
						std::vector<int>split_to_pos_list;
						if (to_index == 0)
						{
							split_to_pos_list.clear();
							for (int tmp_to_pos = to_up_bound[to_index] - 1; tmp_to_pos >= to_low_bound[to_index]; tmp_to_pos--)
							{
								if (count == 0)
								{
									if (!(left_clipped_reads_dict.find(ref_2) == left_clipped_reads_dict.end()) && !(left_clipped_reads_dict[ref_2].find(tmp_to_pos) == left_clipped_reads_dict[ref_2].end()))
										split_to_pos_list.insert(split_to_pos_list.begin(), tmp_to_pos);
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
								else
								{
									if (!(left_clipped_reads_dict.find(ref_2) == left_clipped_reads_dict.end()) && !(left_clipped_reads_dict[ref_2].find(tmp_to_pos) == left_clipped_reads_dict[ref_2].end()))
										continue;
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
							}
						}
						else
						{
							split_to_pos_list.clear();
							for (int tmp_to_pos = to_low_bound[to_index]; tmp_to_pos < to_up_bound[to_index]; tmp_to_pos++)
							{
								if (count == 0)
								{
									if (!(right_clipped_reads_dict.find(ref_2) == right_clipped_reads_dict.end()) && !(right_clipped_reads_dict[ref_2].find(tmp_to_pos) == right_clipped_reads_dict[ref_2].end()))
										split_to_pos_list.insert(split_to_pos_list.begin(), tmp_to_pos);
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
								else
								{
									if (!(right_clipped_reads_dict.find(ref_2) == right_clipped_reads_dict.end()) && !(right_clipped_reads_dict[ref_2].find(tmp_to_pos) == right_clipped_reads_dict[ref_2].end()))
										continue;
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
							}
						}
						for (int tmp_from_pos = from_up_bound[from_index] - 1; tmp_from_pos >= from_low_bound[from_index]; tmp_from_pos--)
						{
							if (!(left_clipped_reads_dict[ref_1].find(tmp_from_pos) == left_clipped_reads_dict[ref_1].end()))
							{
								auto curTime = std::chrono::system_clock::now();
								std::chrono::duration<double> elapsed_seconds = curTime - startTime;
								if (elapsed_seconds.count() > time_limit && !found_flag)
								{
									break;
								}
								if (found_flag)
									break;
								if (to_index == 0)
								{
									for (int i = 0; static_cast<unsigned int>(i) < split_to_pos_list.size(); i++)
									{
										auto curTime = std::chrono::system_clock::now();
										std::chrono::duration<double> elapsed_seconds = curTime - startTime;
										if (elapsed_seconds.count() > time_limit && !found_flag)
										{
											break;
										}
										if (found_flag)
											break;
										int tmp_to_pos = split_to_pos_list[i];
										for (int j = 0; static_cast<unsigned int>(j) < left_clipped_reads_dict[ref_1][tmp_from_pos].size(); j++)
										{
											std::string s1 = left_clipped_reads_dict[ref_1][tmp_from_pos][j];

											int up_bound = tmp_to_pos + static_cast<int>(s1.size());
											if (up_bound > RefSeqLength[ref_2])
												up_bound = RefSeqLength[ref_2];
											if (tmp_to_pos >= up_bound)
												continue;
											std::string tmp_s2 = rg->QueryRegion(ref_2, tmp_to_pos, up_bound - 1);
											std::string s2 = getMateSeq(tmp_s2);
											double similar_degree12 = 0;
											if (s1.size() < 5 || s2.size() < 5)
												similar_degree12 = 0;
											else
											{
												int matches = sswAlignment(s2, s1);
												similar_degree12 = (double)matches / s1.size();
											}
											if (!(left_clipped_reads_dict.find(ref_2) == left_clipped_reads_dict.end()) && !(left_clipped_reads_dict[ref_2].find(tmp_to_pos) == left_clipped_reads_dict[ref_2].end()))
											{
												for (int k = 0; static_cast<unsigned int>(k) < left_clipped_reads_dict[ref_2][tmp_to_pos].size(); k++)
												{
													std::string s4 = left_clipped_reads_dict[ref_2][tmp_to_pos][k];
													if (similar_degree12 < 0.8)
													{
														if (s4.size() <= max_match_split_length)
															continue;
													}
													if ((s4.size() + s2.size() <= max_match_split_length) || s4.size() == 0)
														continue;
													int up_bound = tmp_from_pos + static_cast<int>(s4.size());
													if (up_bound > RefSeqLength[ref_1])
														up_bound = RefSeqLength[ref_1];
													if (tmp_from_pos >= up_bound)
														continue;
													std::string tmp_s3 = rg->QueryRegion(ref_1, tmp_from_pos, up_bound - 1);
													std::string s3 = getMateSeq(tmp_s3);
													double similar_degree34 = 0;
													if (s3.size() < 5 || s4.size() < 5)
														similar_degree34 = 0;
													else
													{
														int matches = sswAlignment(s3, s4);
														similar_degree34 = (double)matches / s4.size();
													}
													if (similar_degree12 >= 0.8 && similar_degree34 < 0.8 && s2.size() > max_match_split_length)
													{
														max_match_split_length = s2.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("left");
															to_side = std::string("left");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s1;
															ref_str = s2;
															similarity = similar_degree12;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 < 0.8 && similar_degree34 >= 0.8 && s4.size() > max_match_split_length)
													{
														max_match_split_length = s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("left");
															to_side = std::string("left");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s4;
															ref_str = s3;
															similarity = similar_degree34;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 >= 0.8 && similar_degree34 >= 0.8 && s2.size() + s4.size() > max_match_split_length)
													{
														max_match_split_length = s2.size() + s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("left");
															to_side = std::string("left");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s4;
															ref_str = s3;
															similarity = similar_degree34;
															found_flag = true;
															break;
														}
													}
												}
											}
											else
											{
												if (similar_degree12 >= 0.8 && s2.size() > max_match_split_length)
												{
													max_match_split_length = s2.size();
													if (max_match_split_length > threshold)
													{
														from_side = std::string("left");
														to_side = std::string("left");
														to_pos = tmp_to_pos;
														from_pos = tmp_from_pos;
														read_str = s1;
														ref_str = s2;
														similarity = similar_degree12;
														found_flag = true;
														break;
													}
												}
											}
										}
									}
								}
								else
								{
									for (int i = 0; static_cast<unsigned int>(i) < split_to_pos_list.size(); i++)
									{
										auto curTime = std::chrono::system_clock::now();
										std::chrono::duration<double> elapsed_seconds = curTime - startTime;
										if (elapsed_seconds.count() > time_limit && !found_flag)
										{
											break;
										}
										if (found_flag)
											break;
										int tmp_to_pos = split_to_pos_list[i];
										for (int j = 0; static_cast<unsigned int>(j) < left_clipped_reads_dict[ref_1][tmp_from_pos].size(); j++)
										{
											std::string s1 = left_clipped_reads_dict[ref_1][tmp_from_pos][j];
											int tmp_low_bound = tmp_to_pos - static_cast<int>(s1.size());
											if (tmp_low_bound < 0)
												tmp_low_bound = 0;
											if (tmp_low_bound >= tmp_to_pos)
												continue;
											std::string s2 = rg->QueryRegion(ref_2, tmp_low_bound, tmp_to_pos - 1);
											double similar_degree12 = 0;
											if (s1.size() < 5 || s2.size() < 5)
												similar_degree12 = 0;
											else
											{
												int matches = sswAlignment(s2, s1);
												similar_degree12 = (double)matches / s1.size();
											}
											if (!(right_clipped_reads_dict.find(ref_2) == right_clipped_reads_dict.end()) && !(right_clipped_reads_dict[ref_2].find(tmp_to_pos) == right_clipped_reads_dict[ref_2].end()))
											{
												for (int k = 0; static_cast<unsigned int>(k) < right_clipped_reads_dict[ref_2][tmp_to_pos].size(); k++)
												{
													std::string s4 = right_clipped_reads_dict[ref_2][tmp_to_pos][k];
													if (similar_degree12 < 0.8)
													{
														if (s4.size() <= max_match_split_length)
															continue;
													}
													if ((s2.size() + s4.size() <= max_match_split_length) || s4.size() == 0)
														continue;
													int up_bound = tmp_from_pos + static_cast<int>(s4.size());
													if (up_bound > RefSeqLength[ref_1])
														up_bound = RefSeqLength[ref_1];
													if (tmp_from_pos >= up_bound)
														continue;
													std::string s3 = rg->QueryRegion(ref_1, tmp_from_pos, up_bound - 1);
													double similar_degree34 = 0;
													if (s3.size() < 5 || s4.size() < 5)
														similar_degree34 = 0;
													else
													{
														int matches = sswAlignment(s3, s4);
														similar_degree34 = (double)matches / s4.size();
													}
													if (similar_degree12 >= 0.8 && similar_degree34 < 0.8 && s2.size() > max_match_split_length)
													{
														max_match_split_length = s2.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("left");
															to_side = std::string("right");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s1;
															ref_str = s2;
															similarity = similar_degree12;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 < 0.8 && similar_degree34 >= 0.8 && s4.size() > max_match_split_length)
													{
														max_match_split_length = s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("left");
															to_side = std::string("right");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s3;
															ref_str = s4;
															similarity = similar_degree34;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 >= 0.8 && similar_degree34 >= 0.8 && s2.size() + s4.size() > max_match_split_length)
													{
														max_match_split_length = s2.size() + s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("left");
															to_side = std::string("right");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s1;
															ref_str = s2;
															similarity = similar_degree12;
															found_flag = true;
															break;
														}
													}
												}
											}
											else
											{
												if (similar_degree12 >= 0.8 && s2.size() > max_match_split_length)
												{
													max_match_split_length = s2.size();
													if (max_match_split_length > threshold)
													{
														from_side = std::string("left");
														to_side = std::string("right");
														to_pos = tmp_to_pos;
														from_pos = tmp_from_pos;
														read_str = s1;
														ref_str = s2;
														similarity = similar_degree12;
														found_flag = true;
														break;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
				else
				{
					if (!(right_clipped_reads_dict.find(ref_1) == right_clipped_reads_dict.end()))
					{
						auto curTime = std::chrono::system_clock::now();
						std::chrono::duration<double> elapsed_seconds = curTime - startTime;
						if (elapsed_seconds.count() > time_limit && !found_flag)
						{
							break;
						}
						if (found_flag)
							break;
						std::vector<int>split_to_pos_list;
						if (to_index == 0)
						{
							split_to_pos_list.clear();
							for (int tmp_to_pos = to_up_bound[to_index] - 1; tmp_to_pos >= to_low_bound[to_index]; tmp_to_pos--)
							{
								if (count == 0)
								{
									if (!(left_clipped_reads_dict.find(ref_2) == left_clipped_reads_dict.end()) && !(left_clipped_reads_dict[ref_2].find(tmp_to_pos) == left_clipped_reads_dict[ref_2].end()))
										split_to_pos_list.insert(split_to_pos_list.begin(), tmp_to_pos);
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
								else
								{
									if (!(left_clipped_reads_dict.find(ref_2) == left_clipped_reads_dict.end()) && !(left_clipped_reads_dict[ref_2].find(tmp_to_pos) == left_clipped_reads_dict[ref_2].end()))
										continue;
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
							}
						}
						else
						{
							split_to_pos_list.clear();
							for (int tmp_to_pos = to_low_bound[to_index]; tmp_to_pos < to_up_bound[to_index]; tmp_to_pos++)
							{
								if (count == 0)
								{
									if (!(right_clipped_reads_dict.find(ref_2) == right_clipped_reads_dict.end()) && !(right_clipped_reads_dict[ref_2].find(tmp_to_pos) == right_clipped_reads_dict[ref_2].end()))
										split_to_pos_list.insert(split_to_pos_list.begin(), tmp_to_pos);
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
								else
								{
									if (!(right_clipped_reads_dict.find(ref_2) == right_clipped_reads_dict.end()) && !(right_clipped_reads_dict[ref_2].find(tmp_to_pos) == right_clipped_reads_dict[ref_2].end()))
										continue;
									else
										split_to_pos_list.push_back(tmp_to_pos);
								}
							}
						}
						for (int tmp_from_pos = from_low_bound[from_index]; tmp_from_pos < from_up_bound[from_index]; tmp_from_pos++)
						{
							if (!(right_clipped_reads_dict[ref_1].find(tmp_from_pos) == right_clipped_reads_dict[ref_1].end()))
							{
								auto curTime = std::chrono::system_clock::now();
								std::chrono::duration<double> elapsed_seconds = curTime - startTime;
								if (elapsed_seconds.count() > time_limit && !found_flag)
								{
									break;
								}
								if (found_flag)
									break;
								if (to_index == 0)
								{
									for (int i = 0; static_cast<unsigned int>(i) < split_to_pos_list.size(); i++)
									{
										auto curTime = std::chrono::system_clock::now();
										std::chrono::duration<double> elapsed_seconds = curTime - startTime;
										if (elapsed_seconds.count() > time_limit && !found_flag)
										{
											break;
										}
										if (found_flag)
											break;
										int tmp_to_pos = split_to_pos_list[i];
										for (int j = 0; static_cast<unsigned int>(j) < right_clipped_reads_dict[ref_1][tmp_from_pos].size(); j++)
										{
											std::string s1 = right_clipped_reads_dict[ref_1][tmp_from_pos][j];
											int up_bound = tmp_to_pos + static_cast<int>(s1.size());
											if (up_bound > RefSeqLength[ref_2])
												up_bound = RefSeqLength[ref_2];
											if (tmp_to_pos >= up_bound)
												continue;
											std::string s2 = rg->QueryRegion(ref_2, tmp_to_pos, up_bound - 1);
											double similar_degree12 = 0;
											if (s1.size() < 5 || s2.size() < 5)
												similar_degree12 = 0;
											else
											{
												int matches = sswAlignment(s2, s1);
												similar_degree12 = (double)matches / s1.size();
											}
											if (!(left_clipped_reads_dict.find(ref_2) == left_clipped_reads_dict.end()) && !(left_clipped_reads_dict[ref_2].find(tmp_to_pos) == left_clipped_reads_dict[ref_2].end()))
											{
												for (int k = 0; static_cast<unsigned int>(k) < left_clipped_reads_dict[ref_2][tmp_to_pos].size(); k++)
												{
													std::string s4 = left_clipped_reads_dict[ref_2][tmp_to_pos][k];
													int sub_low_bound = tmp_from_pos - s4.size();
													if (sub_low_bound < 0)
														sub_low_bound = 0;
													if (similar_degree12 < 0.8)
													{
														if (tmp_from_pos - sub_low_bound < static_cast<int>(max_match_split_length))
															continue;
													}
													if ((tmp_from_pos - sub_low_bound + s2.size() <= max_match_split_length) || sub_low_bound == tmp_from_pos)
														continue;
													if (sub_low_bound >= tmp_from_pos)
														continue;
													std::string s3 = rg->QueryRegion(ref_1, sub_low_bound, tmp_from_pos - 1);
													double similar_degree34 = 0;
													if (s3.size() < 5 || s4.size() < 5)
														similar_degree34 = 0;
													else
													{
														int matches = sswAlignment(s3, s4);
														similar_degree34 = (double)matches / s4.size();
													}
													if (similar_degree12 >= 0.8 && similar_degree34 < 0.8 && s2.size() > max_match_split_length)
													{
														max_match_split_length = s2.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("right");
															to_side = std::string("left");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s1;
															ref_str = s2;
															similarity = similar_degree12;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 < 0.8 && similar_degree34 >= 0.8 && s4.size() > max_match_split_length)
													{
														max_match_split_length = s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("right");
															to_side = std::string("left");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s3;
															ref_str = s4;
															similarity = similar_degree34;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 >= 0.8 && similar_degree34 >= 0.8 && s2.size() + s4.size() > max_match_split_length)
													{
														max_match_split_length = s2.size() + s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("right");
															to_side = std::string("left");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s1;
															ref_str = s2;
															similarity = similar_degree12;
															found_flag = true;
															break;
														}
													}
												}
											}
											else
											{
												if (similar_degree12 >= 0.8 && s2.size() > max_match_split_length)
												{
													max_match_split_length = s2.size();
													if (max_match_split_length > threshold)
													{
														from_side = std::string("right");
														to_side = std::string("left");
														to_pos = tmp_to_pos;
														from_pos = tmp_from_pos;
														read_str = s1;
														ref_str = s2;
														similarity = similar_degree12;
														found_flag = true;
														break;
													}
												}
											}
										}
									}
								}
								else
								{
									for (int i = 0; static_cast<unsigned int>(i) < split_to_pos_list.size(); i++)
									{
										auto curTime = std::chrono::system_clock::now();
										std::chrono::duration<double> elapsed_seconds = curTime - startTime;
										if (elapsed_seconds.count() > time_limit && !found_flag)
										{
											break;
										}
										if (found_flag)
											break;
										int tmp_to_pos = split_to_pos_list[i];
										for (int j = 0; static_cast<unsigned int>(j) < right_clipped_reads_dict[ref_1][tmp_from_pos].size(); j++)
										{
											std::string s1 = right_clipped_reads_dict[ref_1][tmp_from_pos][j];
											int tmp_low_bound = tmp_to_pos - static_cast<int>(s1.size());
											if (tmp_low_bound < 0)
												tmp_low_bound = 0;
											if (tmp_low_bound >= tmp_to_pos)
												continue;
											std::string tmp_s2 = rg->QueryRegion(ref_2, tmp_low_bound, tmp_to_pos - 1);
											std::string s2 = getMateSeq(tmp_s2);
											double similar_degree12 = 0;
											if (s1.size() < 5 || s2.size() < 5)
												similar_degree12 = 0;
											else
											{
												int matches = sswAlignment(s2, s1);
												similar_degree12 = (double)matches / s1.size();
											}
											if (!(right_clipped_reads_dict.find(ref_2) == right_clipped_reads_dict.end()) && !(right_clipped_reads_dict[ref_2].find(tmp_to_pos) == right_clipped_reads_dict[ref_2].end()))
											{
												for (int k = 0; static_cast<unsigned int>(k) < right_clipped_reads_dict[ref_2][tmp_to_pos].size(); k++)
												{
													std::string s4 = right_clipped_reads_dict[ref_2][tmp_to_pos][k];
													int sub_low_bound = tmp_from_pos - s4.size();
													if (sub_low_bound < 0)
														sub_low_bound = 0;
													if (similar_degree12 < 0.8)
													{
														if (tmp_from_pos - sub_low_bound <= static_cast<int>(max_match_split_length))
															continue;
													}
													if (tmp_from_pos - sub_low_bound + s2.size() <= max_match_split_length)
														continue;
													if (sub_low_bound >= tmp_from_pos)
														continue;
													std::string tmp_s3 = rg->QueryRegion(ref_1, sub_low_bound, tmp_from_pos - 1);
													std::string s3 = getMateSeq(tmp_s3);
													double similar_degree34 = 0;
													if (s3.size() < 5 || s4.size() < 5)
														similar_degree34 = 0;
													else
													{
														int matches = sswAlignment(s3, s4);
														similar_degree34 = (double)matches / s4.size();
													}
													if (similar_degree12 >= 0.8 && similar_degree34 < 0.8 && s2.size() > max_match_split_length)
													{
														max_match_split_length = s2.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("right");
															to_side = std::string("right");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s1;
															ref_str = s2;
															similarity = similar_degree12;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 < 0.8 && similar_degree34 >= 0.8 && s4.size() > max_match_split_length)
													{
														max_match_split_length = s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("right");
															to_side = std::string("right");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s3;
															ref_str = s4;
															similarity = similar_degree34;
															found_flag = true;
															break;
														}
													}
													if (similar_degree12 >= 0.8 && similar_degree34 >= 0.8 && s2.size() + s4.size() > max_match_split_length)
													{
														max_match_split_length = s2.size() + s4.size();
														if (max_match_split_length > threshold)
														{
															from_side = std::string("right");
															to_side = std::string("right");
															to_pos = tmp_to_pos;
															from_pos = tmp_from_pos;
															read_str = s1;
															ref_str = s2;
															similarity = similar_degree12;
															found_flag = true;
															break;
														}
													}
												}
											}
											else
											{
												if (similar_degree12 >= 0.8 && s2.size() > max_match_split_length)
												{
													max_match_split_length = s2.size();
													if (max_match_split_length > threshold)
													{
														from_side = std::string("right");
														to_side = std::string("right");
														to_pos = tmp_to_pos;
														from_pos = tmp_from_pos;
														read_str = s1;
														ref_str = s2;
														similarity = similar_degree12;
														found_flag = true;
														break;
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
			if (max_match_split_length > threshold)
			{
				accBkp tmp_bkp;
				auto curTime = std::chrono::system_clock::now();
				std::chrono::duration<double> elapsed_seconds = curTime - startTime;
				std::string interval = std::to_string(elapsed_seconds.count());
				tmp_bkp.interval = interval;
				tmp_bkp.from_ref = ref_1;
				tmp_bkp.to_ref = ref_2;
				tmp_bkp.from_bkp = std::to_string(from_pos);
				tmp_bkp.to_bkp = std::to_string(to_pos);
				tmp_bkp.from_side = from_side;
				tmp_bkp.to_side = to_side;
				tmp_bkp.read_str = read_str;
				tmp_bkp.ref_str = ref_str;
				tmp_bkp.similarity = std::to_string(similarity);
				if (reverse)
					tmp_bkp.reverse = "true";
				else
					tmp_bkp.reverse = "false";
				better_breakpoint_list.push_back(tmp_bkp);
				break;
			}
		}
	}
	
	g_mutex.lock();
	delete rg;

	std::ofstream out(output, std::ofstream::app);
	for (int i = 0; static_cast<unsigned int>(i) < better_breakpoint_list.size(); i++)
	{
		accBkp bkp = better_breakpoint_list[i];

		std::string out_bkp = bkp.from_ref + ',' + bkp.from_bkp + ',' + bkp.to_ref + ',' + bkp.to_bkp + ',' + bkp.from_side + ',' + bkp.to_side + ',' + bkp.reverse + ','+ bkp.read_str+','+ bkp.ref_str+','+ bkp.similarity+','+ bkp.interval+'\n';
		out << out_bkp;
	}
	out.close();
	g_mutex.unlock();
}



void GetAccBkp::spawnThreads(int n, std::vector<rawBkp>raw_breakpoints_list, std::unordered_map < std::string, std::unordered_map<int, std::vector < std::string >> > left_clipped_reads_dict,
	std::unordered_map<std::string, std::unordered_map<int, std::vector < std::string >>>right_clipped_reads_dict, std::string output, unsigned int threshold, std::string ref_file, int rlen, int insert_size)
{
	if (raw_breakpoints_list.size() <200)
	{
		n = 2;
	}
	int num_bkp_per = (int)(raw_breakpoints_list.size() / n) + 1;
	std::vector<std::thread>threads;
	for (int i = 0; i < n; i++)
	{
		int start = i * num_bkp_per;
		int end;
		if (i != n - 1)
			end = (i + 1)*num_bkp_per;
		else
			end = raw_breakpoints_list.size();
		threads.push_back(std::thread(getAccurateBreakpoints, start, end, std::ref(raw_breakpoints_list), std::ref(left_clipped_reads_dict), std::ref(right_clipped_reads_dict), output, threshold, ref_file, rlen, insert_size));
	}
	for (auto&th : threads)
		th.join();

}