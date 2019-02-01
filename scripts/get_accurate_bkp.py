import sys
sys.path.append('/home/lichen/software/python36/lib/python3.6/site-packages/parasail/')
import os
import re
import parasail
import pysam
import multiprocessing
import argparse
import math


def readFilter(read):
    return (read.is_proper_pair and
            read.is_paired and
            read.tlen > 0 and
            read.tlen < 1000 and
            not read.is_supplementary and
            not read.is_duplicate and
            not read.is_unmapped and
            not read.mate_is_unmapped)


def getInsertSize(unique_bamfile):
    read_length_list = []
    insert_size_list = []
    for read in unique_bamfile:
        if readFilter(read):
            insert_size_list.append(read.tlen)
            read_length_list.append(len(read.query_sequence))
    read_length = int(sum(read_length_list) / len(read_length_list))
    mean = float(sum(insert_size_list)) / len(insert_size_list)
    sdev = math.sqrt(float(sum([(x - mean)**2 for x in insert_size_list])) / (len(insert_size_list) - 1))
    return mean, sdev, read_length


def getSplitReads(split_bamfile):
    split_reads   =   {}
    for read in split_bamfile: 
        if read.is_unmapped ==  False:
            if read.qname not in split_reads:
                ls = []
                ls.append(read)
                dict_tmp    =   {read.qname:ls}
                split_reads.update(dict_tmp)
            else:
                split_reads.get(read.qname).append(read)
    return split_reads

def indexSplitReadsOnRef(split_reads):
    ref_split_reads = {}
    for key in split_reads:
        if len(split_reads.get(key)) > 0:
            for i in range(0,len(split_reads.get(key))):
                if split_reads.get(key)[i].reference_name not in ref_split_reads:
                    ls = []
                    ls.append(split_reads.get(key)[i])
                    dict_tmp = {split_reads.get(key)[i].reference_name:ls}
                    ref_split_reads.update(dict_tmp)
                else:
                    ref_split_reads.get(split_reads.get(key)[i].reference_name).append(split_reads.get(key)[i])
    return ref_split_reads


def indexSplitReadOnPos(ref_split_reads):
    ref_pos_split_reads_zero = {}
    ref_pos_split_reads_nonzero = {}
    for key in ref_split_reads:
        tmp_dict_zero = {}
        tmp_dict_nonzero = {}
        for i in range(0, len(ref_split_reads.get(key))):
            read = ref_split_reads.get(key)[i]
            if len(read.query_sequence[read.query_alignment_end:]) != 0:
                bkp_pos = read.reference_start + read.query_alignment_length
                if bkp_pos not in tmp_dict_zero:
                    ls = []
                    ls.append(read)
                    buf_dict = {bkp_pos : ls}
                    tmp_dict_zero.update(buf_dict)
                else:
                    tmp_dict_zero.get(bkp_pos).append(read)
            if len(read.query_sequence[:read.query_alignment_start]) != 0:
                bkp_pos = read.reference_start
                if bkp_pos not in tmp_dict_nonzero:
                    ls = []
                    ls.append(read)
                    buf_dict = {bkp_pos : ls}
                    tmp_dict_nonzero.update(buf_dict)
                else:
                    tmp_dict_nonzero.get(bkp_pos).append(read)
        if len(tmp_dict_zero) > 0:
            tmp_ref_dict_zero = {key : tmp_dict_zero}
            ref_pos_split_reads_zero.update(tmp_ref_dict_zero)
        if len(tmp_dict_nonzero) > 0:
            tmp_ref_dict_nonzero = {key : tmp_dict_nonzero}
            ref_pos_split_reads_nonzero.update(tmp_ref_dict_nonzero)
    return ref_pos_split_reads_zero, ref_pos_split_reads_nonzero



def get_mate_seq(sequence):
    mate_seq = ''
    for i in range(len(sequence)-1, -1, -1):
        if sequence[i] == 'A':
            mate_seq += 'T'
        elif sequence[i] == 'T':
            mate_seq += 'A'
        elif sequence[i] == 'G':
            mate_seq += 'C'
        elif sequence[i] == 'C':
            mate_seq += 'G'
    return mate_seq


def worker(candidate_breakpoints, ref_pos_split_reads_zero, ref_pos_split_reads_nonzero, out_put_file, threshold, ref_genome):
    better_breakpoint = []
    for cbk in candidate_breakpoints:
        key = cbk[0]
        sub_key = cbk[1]
        from_seg = cbk[-1]
        if from_seg[2] == key:
            continue
        if from_seg[0] < from_seg[1]:
            left_low = from_seg[0] - insert_size
            right_up = from_seg[1] + rlen + insert_size
            left_up = from_seg[1] + threshold
            right_low = from_seg[1] - threshold
        else:
            left_low = from_seg[1] - insert_size
            right_up = from_seg[0] + rlen + insert_size
            left_up = from_seg[0] + threshold
            right_low = from_seg[0] - threshold
        if left_low < 1:
            left_low = 1
        if right_low < 0:
            right_low = 1
        if left_low > left_up or right_low > right_up:
            continue
        low_bound = [left_low, right_low]
        up_bound = [left_up, right_up]
        query_key = from_seg[2]
        reverse = from_seg[-1]
        if from_seg[4] < from_seg[5]:
            to_left_low = from_seg[4] - insert_size
            to_right_up = from_seg[5] + rlen + insert_size
            to_left_up = from_seg[5] + threshold
            to_right_low = from_seg[5] - threshold
        else:
            to_left_low = from_seg[5] - insert_size
            to_right_up = from_seg[4] + rlen + insert_size
            to_left_up = from_seg[4] + threshold
            to_right_low = from_seg[4] - threshold
        if to_left_low < 0:
            to_left_low = 1
        if to_right_low < 0:
            to_right_low = 1
        if to_left_low > to_left_up or to_right_low > to_right_up:
            continue
        to_low_bound = [to_left_low, to_right_low]
        to_up_bound = [to_left_up, to_right_up]
        to_pos = from_seg[3]
        from_pos = sub_key
        max_match_split_length = 0
        from_side = 0
        to_side = 0
        for to_index in range(0, 2):
            if reverse == False:
                if to_index == 0:
                    from_index = 1
                else:
                    from_index = 0
            else:
                if to_index == 0:
                    from_index = 0
                else:
                    from_index = 1
            for tmp_to_pos in range(to_low_bound[to_index], to_up_bound[to_index]):
                for tmp_from_pos in range(low_bound[from_index], up_bound[from_index]):
                    if from_index == 0:
                        if to_index == 0:
                            if key in ref_pos_split_reads_nonzero and tmp_from_pos in ref_pos_split_reads_nonzero[key]:
                                for j in range(len(ref_pos_split_reads_nonzero[key][tmp_from_pos])):
                                    soft_clip_read = ref_pos_split_reads_nonzero[key][tmp_from_pos][j]
                                    s1 = soft_clip_read.query_sequence[:soft_clip_read.query_alignment_start]
                                    s2 = ref_genome.fetch(query_key, tmp_to_pos, tmp_to_pos + len(s1))
                                    s2 = get_mate_seq(s2)
                                    al = parasail.sw_stats_striped_64(s1, s2, 2, -2, parasail.pam100)
                                    similar_degree12 = al.matches / len(s1)
                                    if query_key in ref_pos_split_reads_nonzero and tmp_to_pos in ref_pos_split_reads_nonzero[query_key]:
                                        for soft_clip_read in ref_pos_split_reads_nonzero.get(query_key).get(tmp_to_pos):
                                            s4 = soft_clip_read.query_sequence[:soft_clip_read.query_alignment_start]
                                            if similar_degree12 < 0.8:
                                                if len(s4) <= max_match_split_length:
                                                    continue
                                            if len(s4) + len(s2) <= max_match_split_length:
                                                continue
                                            s3 = ref_genome.fetch(key, tmp_from_pos, tmp_from_pos + len(s4))
                                            s3 = get_mate_seq(s3)
                                            al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                            if len(s4) == 0:
                                                continue
                                            similar_degree34 = al.matches / len(s4)
                                            if similar_degree12 >= 0.8 and similar_degree34 < 0.8 and len(s2) > max_match_split_length:
                                                max_match_split_length = len(s2)
                                                to_side = 'left'
                                                from_side = 'left'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 < 0.8 and similar_degree34 >= 0.8 and len(s4) > max_match_split_length:
                                                max_match_split_length = len(s4)
                                                to_side = 'left'
                                                from_side = 'left'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 >= 0.8 and similar_degree34 >= 0.8 and len(s2) + len(s4) > max_match_split_length:
                                                max_match_split_length = len(s2) + len(s4)
                                                to_side = 'left'
                                                from_side = 'left'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                    else:
                                        if similar_degree12 >= 0.8 and len(s2) > max_match_split_length:
                                            to_side = 'left'
                                            from_side = 'left'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s2)
                            else:
                                if query_key in ref_pos_split_reads_nonzero and tmp_to_pos in ref_pos_split_reads_nonzero[query_key]:
                                    for soft_clip_read in ref_pos_split_reads_nonzero.get(query_key).get(tmp_to_pos):
                                        s4 = soft_clip_read.query_sequence[:soft_clip_read.query_alignment_start]
                                        if len(s4) <= max_match_split_length:
                                            continue
                                        s3 = ref_genome.fetch(key, tmp_from_pos, tmp_from_pos + len(s4))
                                        s3 = get_mate_seq(s3)
                                        al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                        if len(s4) == 0:
                                            continue
                                        similar_degree34 = al.matches / len(s4)
                                        if similar_degree34 >= 0.8:
                                            to_side = 'left'
                                            from_side = 'left'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s4)
                        else:
                            if key in ref_pos_split_reads_nonzero and tmp_from_pos in ref_pos_split_reads_nonzero[key]:
                                for j in range(0, len(ref_pos_split_reads_nonzero.get(key).get(tmp_from_pos))):
                                    soft_clip_read = ref_pos_split_reads_nonzero.get(key).get(tmp_from_pos)[j]
                                    s1 = soft_clip_read.query_sequence[:soft_clip_read.query_alignment_start]
                                    tmp_low_bound = tmp_to_pos - len(s1)
                                    if tmp_low_bound < 0:
                                        tmp_low_bound = 0
                                    if tmp_low_bound > tmp_to_pos:
                                        continue
                                    s2 = ref_genome.fetch(query_key, tmp_low_bound, tmp_to_pos)
                                    al = parasail.sw_stats_striped_64(s1, s2, 2, -2, parasail.pam100)
                                    if len(s2) == 0:
                                        continue
                                    similar_degree12 = al.matches / len(s2)
                                    if query_key in ref_pos_split_reads_zero and tmp_to_pos in ref_pos_split_reads_zero[query_key]:
                                        for soft_clip_read in ref_pos_split_reads_zero.get(query_key).get(tmp_to_pos):
                                            s4 = soft_clip_read.query_sequence[soft_clip_read.query_alignment_end:]
                                            if similar_degree12 < 0.8:
                                                if len(s4) <= max_match_split_length:
                                                    continue
                                            if len(s4) + len(s2) <= max_match_split_length:
                                                continue
                                            s3 = ref_genome.fetch(key, tmp_from_pos, tmp_from_pos + len(s4))
                                            al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                            if len(s4) == 0:
                                                continue
                                            similar_degree34 = al.matches / len(s4)
                                            if similar_degree12 >= 0.8 and similar_degree34 < 0.8 and len(s2) > max_match_split_length:
                                                max_match_split_length = len(s2)
                                                to_side = 'right'
                                                from_side = 'left'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 < 0.8 and similar_degree34 >= 0.8 and len(s4) > max_match_split_length:
                                                max_match_split_length = len(s4)
                                                to_side = 'right'
                                                from_side = 'left'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 >= 0.8 and similar_degree34 >= 0.8 and len(s2) + len(s4) > max_match_split_length:
                                                max_match_split_length = len(s2) + len(s4)
                                                to_side = 'right'
                                                from_side = 'left'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                    else:
                                        if similar_degree12 >= 0.8 and len(s2) > max_match_split_length:
                                            to_side = 'right'
                                            from_side = 'left'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s2)
                            else:
                                if query_key in ref_pos_split_reads_zero and tmp_to_pos in ref_pos_split_reads_zero[query_key]:
                                    for soft_clip_read in ref_pos_split_reads_zero.get(query_key).get(tmp_to_pos):
                                        s4 = soft_clip_read.query_sequence[soft_clip_read.query_alignment_end:]
                                        if len(s4) <= max_match_split_length:
                                            continue
                                        s3 = ref_genome.fetch(key, tmp_from_pos, tmp_from_pos + len(s4))
                                        al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                        if len(s4) == 0:
                                            continue
                                        similar_degree34 = al.matches / len(s4)
                                        if similar_degree34 >= 0.8:
                                            to_side = 'right'
                                            from_side = 'left'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s4)
                    else:
                        if to_index == 0:
                            if key in ref_pos_split_reads_zero and tmp_from_pos in ref_pos_split_reads_zero[key]:
                                for soft_clip_read in ref_pos_split_reads_zero.get(key).get(tmp_from_pos):
                                    s1 = soft_clip_read.query_sequence[soft_clip_read.query_alignment_end:]
                                    s2 = ref_genome.fetch(query_key, tmp_to_pos, tmp_to_pos + len(s1))
                                    al = parasail.sw_stats_striped_64(s1, s2, 2, -2, parasail.pam100)
                                    if len(s2) == 0:
                                        continue
                                    similar_degree12 = al.matches / len(s2)
                                    if query_key in ref_pos_split_reads_nonzero and tmp_to_pos in ref_pos_split_reads_nonzero[query_key]:
                                        for soft_clip_read in ref_pos_split_reads_nonzero.get(query_key).get(tmp_to_pos):
                                            s4 = soft_clip_read.query_sequence[:soft_clip_read.query_alignment_start]
                                            sub_low_bound = tmp_from_pos - len(s4)
                                            if sub_low_bound < 0:
                                                sub_low_bound = 0
                                            if similar_degree12 < 0.8:
                                                if tmp_from_pos - sub_low_bound < max_match_split_length:
                                                    continue
                                            if tmp_from_pos - sub_low_bound + len(s2) <= max_match_split_length:
                                                continue
                                            s3 = ref_genome.fetch(key, sub_low_bound, tmp_from_pos)
                                            al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                            if len(s3) == 0:
                                                continue
                                            similar_degree34 = al.matches / len(s3)
                                            if similar_degree12 >= 0.8 and similar_degree34 < 0.8 and len(s2) > max_match_split_length:
                                                max_match_split_length = len(s2)
                                                to_side = 'left'
                                                from_side = 'right'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 < 0.8 and similar_degree34 >= 0.8 and len(s4) > max_match_split_length:
                                                max_match_split_length = len(s4)
                                                to_side = 'left'
                                                from_side = 'right'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 >= 0.8 and similar_degree34 >= 0.8 and len(s2) + len(s4) > max_match_split_length:
                                                max_match_split_length = len(s2) + len(s4)
                                                to_side = 'left'
                                                from_side = 'right'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                    else:
                                        if similar_degree12 >= 0.8 and len(s2) > max_match_split_length:
                                            to_side = 'left'
                                            from_side = 'right'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s2)
                            else:
                                if query_key in ref_pos_split_reads_nonzero and tmp_to_pos in ref_pos_split_reads_nonzero[query_key]:
                                    for soft_clip_read in ref_pos_split_reads_nonzero.get(query_key).get(tmp_to_pos):
                                        s4 = soft_clip_read.query_sequence[:soft_clip_read.query_alignment_start]
                                        sub_low_bound = tmp_from_pos - len(s4)
                                        if sub_low_bound < 0:
                                            sub_low_bound = 0
                                        if tmp_from_pos - sub_low_bound < max_match_split_length:
                                            continue
                                        s3 = ref_genome.fetch(key, sub_low_bound, tmp_from_pos)
                                        al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                        if len(s3) == 0:
                                            continue
                                        similar_degree34 = al.matches / len(s3)
                                        if similar_degree34 >= 0.8:
                                            to_side = 'left'
                                            from_side = 'right'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s3)
                        else:
                            if key in ref_pos_split_reads_zero and tmp_from_pos in ref_pos_split_reads_zero[key]:
                                for soft_clip_read in ref_pos_split_reads_zero.get(key).get(tmp_from_pos):
                                    s1 = soft_clip_read.query_sequence[soft_clip_read.query_alignment_end:]
                                    tmp_low_bound = tmp_to_pos - len(s1)
                                    if tmp_low_bound < 0:
                                        tmp_low_bound = 0
                                    if tmp_low_bound > tmp_to_pos:
                                        continue
                                    s2 = ref_genome.fetch(query_key, tmp_low_bound, tmp_to_pos)
                                    if len(s2) == 0:
                                        continue
                                    s2 = get_mate_seq(s2)
                                    al = parasail.sw_stats_striped_64(s1, s2, 2, -2, parasail.pam100)
                                    if len(s2) == 0:
                                        continue
                                    similar_degree12 = al.matches / len(s2)
                                    if query_key in ref_pos_split_reads_zero and tmp_to_pos in ref_pos_split_reads_zero[query_key]:
                                        for sub_soft_clip_read in ref_pos_split_reads_zero.get(query_key).get(tmp_to_pos):
                                            s4 = sub_soft_clip_read.query_sequence[sub_soft_clip_read.query_alignment_end:]
                                            sub_low_bound = tmp_from_pos - len(s4)
                                            if sub_low_bound < 0:
                                                sub_low_bound = 0
                                            if similar_degree12 < 0.8:
                                                if tmp_from_pos - sub_low_bound <= max_match_split_length:
                                                    continue
                                            if tmp_from_pos - sub_low_bound + len(s2) <= max_match_split_length:
                                                continue
                                            s3 = ref_genome.fetch(key, sub_low_bound, tmp_from_pos)
                                            s3 = get_mate_seq(s3)
                                            al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                            if len(s3) == 0:
                                                continue
                                            similar_degree34 = al.matches / len(s3)
                                            if similar_degree12 >= 0.8 and similar_degree34 < 0.8 and len(s2) > max_match_split_length:
                                                max_match_split_length = len(s2)
                                                to_side = 'right'
                                                from_side = 'right'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 < 0.8 and similar_degree34 >= 0.8 and len(s4) > max_match_split_length:
                                                max_match_split_length = len(s4)
                                                to_side = 'right'
                                                from_side = 'right'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                            if similar_degree12 >= 0.8 and similar_degree34 >= 0.8 and len(s2) + len(s4) > max_match_split_length:
                                                max_match_split_length = len(s2) + len(s4)
                                                to_side = 'right'
                                                from_side = 'right'
                                                to_pos = tmp_to_pos
                                                from_pos = tmp_from_pos
                                    else:
                                        if similar_degree12 >= 0.8 and len(s2) > max_match_split_length:
                                            to_side = 'right'
                                            from_side = 'right'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s2)
                            else:
                                if query_key in ref_pos_split_reads_zero and tmp_to_pos in ref_pos_split_reads_zero[query_key]:
                                    for sub_soft_clip_read in ref_pos_split_reads_zero.get(query_key).get(tmp_to_pos):
                                        s4 = sub_soft_clip_read.query_sequence[sub_soft_clip_read.query_alignment_end:]
                                        sub_low_bound = tmp_from_pos - len(s4)
                                        if sub_low_bound < 0:
                                            sub_low_bound = 0
                                        if tmp_from_pos - sub_low_bound <= max_match_split_length:
                                            continue
                                        s3 = ref_genome.fetch(key, sub_low_bound, tmp_from_pos)
                                        s3 = get_mate_seq(s3)
                                        al = parasail.sw_stats_striped_64(s3, s4, 2, -2, parasail.pam100)
                                        if len(s3) == 0:
                                            continue
                                        similar_degree34 = al.matches / len(s3)
                                        if similar_degree34 >= 0.8:
                                            to_side = 'right'
                                            from_side = 'right'
                                            to_pos = tmp_to_pos
                                            from_pos = tmp_from_pos
                                            max_match_split_length = len(s3)
        if max_match_split_length > threshold:
            print(max_match_split_length)
            tmp = [key, from_pos, from_seg[0], from_seg[1], from_seg[2], 
                to_pos, from_seg[4], from_seg[5], from_seg[6], from_side, to_side, str(reverse)]
            better_breakpoint.append(tmp)
    fo = open(out_put_file,  "a")
    for bkp in better_breakpoint:
        key = bkp[0]
        pos1_left = bkp[2]
        pos1_right = bkp[3]
        pos2_left = bkp[6]
        pos2_right = bkp[7]
        to_ref = bkp[4]
        to_pos = bkp[5]
        breakpoint_pos = bkp[1]
        pair_end_coverage = bkp[-4]
        from_side = bkp[-3]
        to_side = bkp[-2]
        reverse = bkp[-1]
        fo.write(key + ',' + str(breakpoint_pos) + ',' + str(pos1_left)+',' + str(pos1_right)+','+ to_ref + ',' + str(to_pos) + ',' + 
            str(pos2_left)+','+str(pos2_right)+ ',' +str(pair_end_coverage) + ',' + from_side + ',' + to_side + ',' + reverse+'\n')
    fo.close()
 


def main():
    split_num = args["t"]
    ref_genome = pysam.Fastafile(args["r"])
    threshold = 15
    out_put_file = args["o"]
    split_bam_name = args["s"]
    split_bamfile = pysam.AlignmentFile(filename = split_bam_name, mode = 'rb')
    split_reads = getSplitReads(split_bamfile)
    ref_split_reads = indexSplitReadsOnRef(split_reads)
    ref_pos_split_reads_zero, ref_pos_split_reads_nonzero = indexSplitReadOnPos(ref_split_reads)
    candidate_breakpoint_file = args["raw_bkp"]
    fi = open(candidate_breakpoint_file, "r")
    clusters_region = {}
    while 1:
        buf = fi.readline()
        buf = buf.strip('\n')
        if not buf:
            break
        buf = re.split(r'[:;,\s]\s*',  buf)
        reverse = False
        if buf[-1] == 'true':
            reverse = True
        if buf[0] not in ref_pos_split_reads_zero and buf[0] not in ref_pos_split_reads_nonzero:
            continue
        if buf[4] not in ref_pos_split_reads_zero and buf[4] not in ref_pos_split_reads_nonzero:
            continue
        if buf[0] not in clusters_region:
            #pos1_left, pos1_right, ref_name2, pos2, pos2_left, pos2_right
            tmp_dict = {int(buf[1]) : {buf[4] : [[int(buf[2]), int(buf[3]), buf[4], int(buf[5]), int(buf[6]), int(buf[7]), int(buf[8]), reverse]]}}
            per_dict = {buf[0] : tmp_dict}
            clusters_region.update(per_dict)
        else:
            found = False
            for candidate_pos in range(int(buf[1]) - int(threshold/2), int(buf[1]) + int(threshold/2)):
                if candidate_pos in clusters_region[buf[0]]:
                    if buf[4] in clusters_region[buf[0]][candidate_pos]:
                        for to_breakpoint in clusters_region[buf[0]][candidate_pos][buf[4]]:
                            if (int(buf[3]) - int(buf[1])) * (to_breakpoint[1] - candidate_pos) >= 0 and abs(int(buf[5]) - to_breakpoint[3]) < threshold and (int(buf[7]) - int(buf[5])) * (to_breakpoint[5] - to_breakpoint[3]) >= 0 and reverse == to_breakpoint[-1]:
                                if (int(buf[2]) - int(buf[1])) * (to_breakpoint[0] - candidate_pos) >= 0 and (int(buf[6]) - int(buf[5])) * (to_breakpoint[4] - to_breakpoint[3]) >= 0:
                                    found = True
                                    break
                        if found is True:
                            break
            if found is False:
                if int(buf[1]) in clusters_region[buf[0]]:
                    if buf[4] not in clusters_region[buf[0]][int(buf[1])]:
                        tmp_dict = {buf[4] : [[int(buf[2]), int(buf[3]), buf[4], int(buf[5]), int(buf[6]), int(buf[7]), int(buf[8]), reverse]]}
                        clusters_region[buf[0]][int(buf[1])].update(tmp_dict)
                    else:
                        clusters_region[buf[0]][int(buf[1])][buf[4]].append([int(buf[2]), int(buf[3]), buf[4], int(buf[5]), int(buf[6]), int(buf[7]), int(buf[8]), reverse])
                else:
                    tmp_dict = {int(buf[1]) : {buf[4] : [[int(buf[2]), int(buf[3]), buf[4], int(buf[5]), int(buf[6]), int(buf[7]), int(buf[8]), reverse]]}}
                    clusters_region[buf[0]].update(tmp_dict)
    clusters_region_list = []
    for key in clusters_region:
        if key not in ref_pos_split_reads_zero and key not in ref_pos_split_reads_nonzero:
            continue
        else:
            for pos_key in clusters_region[key]:
                for to_ref_key in clusters_region[key][pos_key]:
                    for from_seg in clusters_region[key][pos_key][to_ref_key]:
                        tmp = [key, pos_key, from_seg]
                        clusters_region_list.append(tmp)
    num_breakpoint_per = int(len(clusters_region_list) / split_num) + 1
    procs = []
    for i in range(0, split_num):
        if i != split_num - 1:
            candidate_breakpoints = clusters_region_list[i*num_breakpoint_per : (i+1)*num_breakpoint_per]
        else:
            candidate_breakpoints = clusters_region_list[(split_num-1)*num_breakpoint_per:]
        p = multiprocessing.Process(target = worker, args = (candidate_breakpoints, ref_pos_split_reads_zero, ref_pos_split_reads_nonzero, out_put_file, threshold, ref_genome))
        procs.append(p)
        p.start()
        print(p.pid)
    for proc in procs:
        proc.join()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="HGT get accurate breakpoints", add_help=False, usage="%(prog)s [-h] -r genome_dir -id sample_id.txt", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    required = parser.add_argument_group("required arguments")
    optional = parser.add_argument_group("optional arguments")
    required.add_argument("-r", type=str, help="<str> Metagenomic reference", metavar="\b")
    required.add_argument("-o", type=str, help="<str> accurate breakpoints file.", metavar="\b")
    required.add_argument("-s", type=str, help="<str> split reads bam file.", metavar="\b")
    required.add_argument("--raw_bkp", type=str, help="raw breakpoints file.", metavar="\b")
    optional.add_argument("-t", type=int, default=4, help="<int> number of threads", metavar="\b")
    required.add_argument("-u", type=str, help="<str> unique reads bam file.", metavar="\b")
    optional.add_argument("-h", "--help", action="help")
    args = vars(parser.parse_args())
    unique_bam_name = args["u"]
    unique_bamfile = pysam.AlignmentFile(filename = unique_bam_name, mode = 'rb')
    mean, sdev, rlen = getInsertSize(unique_bamfile)
    insert_size = int(mean + sdev)
    rlen = int(rlen)
    print(mean,sdev)
    sys.exit(main())