import sys
import fnmatch
import os
import re
import pysam
import copy
import random
import numpy as np
from scipy import signal
from scipy import ndimage
from itertools import product
import matplotlib.pyplot as plt
import pandas as pd
from gurobipy import *
import argparse



def getGenomeSeqCoverage(file_name, ref_name_list):
    genome_sequence_coverage = {}
    fi = open(file_name, "r")
    while 1:
        buf = fi.readline()
        buf = buf.strip('\n')
        if not buf:
            break
        buf = re.split(r'[:;,\s]\s*', buf)
        ref_name = buf[0]
        if ref_name in ref_name_list:
            left = int(buf[1])
            right = int(buf[2])
            depth = int(buf[3])
            if ref_name not in genome_sequence_coverage:
                per_genome_sequence_coverage = {}
                for i in range(left, right+1):
                    tmp_dict = {i : depth}
                    per_genome_sequence_coverage.update(tmp_dict)
                tmp_dict = {ref_name : per_genome_sequence_coverage}
                genome_sequence_coverage.update(tmp_dict)
            else:
                for i in range(left, right + 1):
                    if i not in genome_sequence_coverage.get(ref_name):
                        tmp_dict = {i : depth}
                        genome_sequence_coverage.get(ref_name).update(tmp_dict)
                    else:
                        new_depth = genome_sequence_coverage[ref_name][i] + depth
                        tmp_dict = {i : new_depth}
                        genome_sequence_coverage[ref_name].update(tmp_dict)
    return genome_sequence_coverage


def modify_edge_breakpoint(filename, genome_sequence_coverage):
    junction_edge_list = []
    change_map_dict = {}
    fi = open(filename, "r")
    while 1:
        buf = fi.readline()
        buf = buf.strip('\n')
        if not buf:
            break
        buf = re.split(r'[:;,\s]\s*', buf)
        from_ref = buf[0]
        from_pos = int(buf[1])
        to_ref = buf[2]
        to_pos = int(buf[3])
        if from_ref not in change_map_dict:
            map_dict = {from_pos : from_pos}
            tmp_dict = {from_ref : map_dict}
            change_map_dict.update(tmp_dict)
        else:
            if from_pos not in change_map_dict[from_ref]:
                found = False
                for pos_key in change_map_dict[from_ref]:
                    if abs(change_map_dict[from_ref][pos_key] - from_pos) < delta:
                        map_dict = {from_pos : change_map_dict[from_ref][pos_key]}
                        change_map_dict[from_ref].update(map_dict)
                        found = True
                        break
                if found is False:
                    map_dict = {from_pos : from_pos}
                    change_map_dict[from_ref].update(map_dict)
        if to_ref not in change_map_dict:
            map_dict = {to_pos : to_pos}
            tmp_dict = {to_ref : map_dict}
            change_map_dict.update(tmp_dict)
        else:
            if to_pos not in change_map_dict[to_ref]:
                found = False
                for pos_key in change_map_dict[to_ref]:
                    if abs(change_map_dict[to_ref][pos_key] - to_pos) < delta:
                        map_dict = {to_pos : change_map_dict[to_ref][pos_key]}
                        change_map_dict[to_ref].update(map_dict)
                        found = True
                        break
                if found is False:
                    map_dict = {to_pos : to_pos}
                    change_map_dict[to_ref].update(map_dict)
    fi.close()
    fi = open(filename, "r")
    while 1:
        buf = fi.readline()
        buf = buf.strip('\n')
        if not buf:
            break
        buf = re.split(r'[:;,\s]\s*', buf)
        from_ref = buf[0]
        from_pos = int(buf[1])
        from_pos = change_map_dict[from_ref][from_pos]
        from_side = buf[4]
        to_ref = buf[2]
        to_pos = int(buf[3])
        to_pos = change_map_dict[to_ref][to_pos]
        to_side = buf[5]
        if buf[6] == 'false':
            reverse = False
        else:
            reverse = True
        if from_ref not in genome_sequence_coverage or to_ref not in genome_sequence_coverage:
            continue
        edge = [from_ref, from_pos, from_side, to_ref, to_pos, to_side, reverse]
        junction_edge_list.append(edge)
    return junction_edge_list


def getDirectedContactBreakpoint(junction_edge_list):
    directed_contact_breakpoint_dict = {}
    for edge in junction_edge_list:
        ref1 = edge[0]
        bkp1 = edge[1]
        bkp1_side = edge[2]
        ref2 = edge[3]
        bkp2 = edge[4]
        bkp2_side = edge[5]
        reverse = edge[-1]
        if bkp1_side == 'right':
            from_ref = ref1
            from_pos = bkp1
            to_ref = ref2
            to_pos = bkp2
            from_side = bkp1_side
            to_side = bkp2_side
        else:
            from_ref = ref2
            from_pos = bkp2
            to_ref = ref1
            to_pos = bkp1
            from_side = bkp2_side
            to_side = bkp1_side
        if from_ref not in directed_contact_breakpoint_dict:
            ls = [[to_pos, from_side, to_side, reverse]]
            sub_ref_dict = {to_ref : ls}
            pos_dict = {from_pos : sub_ref_dict}
            tmp_dict = {from_ref : pos_dict}
            directed_contact_breakpoint_dict.update(tmp_dict)
        else:
            if from_pos not in directed_contact_breakpoint_dict[from_ref]:
                ls = [[to_pos, from_side, to_side, reverse]]
                sub_ref_dict = {to_ref : ls}
                pos_dict = {from_pos : sub_ref_dict}
                directed_contact_breakpoint_dict[from_ref].update(pos_dict)
            else:
                if to_ref not in directed_contact_breakpoint_dict[from_ref][from_pos]:
                    ls = [[to_pos, from_side, to_side, reverse]]
                    sub_ref_dict = {to_ref : ls}
                    directed_contact_breakpoint_dict[from_ref][from_pos].update(sub_ref_dict)
                else:
                    found = False
                    for breakpoint in directed_contact_breakpoint_dict[from_ref][from_pos][to_ref]:
                        if breakpoint[0] == to_pos and breakpoint[-1]  == reverse and breakpoint[1] == from_side and breakpoint[2] == to_side :
                            found = True
                            break
                    if found is False:
                        directed_contact_breakpoint_dict[from_ref][from_pos][to_ref].append([to_pos, from_side, to_side, reverse])
        if reverse is True:
            if bkp1_side == 'right':
                from_ref = ref2
                from_pos = bkp2
                to_ref = ref1
                to_pos = bkp1
                from_side = bkp1_side
                to_side = bkp2_side
            else:
                from_ref = ref1
                from_pos = bkp1
                to_ref = ref2
                to_pos = bkp2
                from_side = bkp1_side
                to_side = bkp2_side
            if from_ref not in directed_contact_breakpoint_dict:
                ls = [[to_pos, from_side, to_side, reverse]]
                sub_ref_dict = {to_ref : ls}
                pos_dict = {from_pos : sub_ref_dict}
                tmp_dict = {from_ref : pos_dict}
                directed_contact_breakpoint_dict.update(tmp_dict)
            else:
                if from_pos not in directed_contact_breakpoint_dict[from_ref]:
                    ls = [[to_pos, from_side, to_side, reverse]]
                    sub_ref_dict = {to_ref : ls}
                    pos_dict = {from_pos : sub_ref_dict}
                    directed_contact_breakpoint_dict[from_ref].update(pos_dict)
                else:
                    if to_ref not in directed_contact_breakpoint_dict[from_ref][from_pos]:
                        ls = [[to_pos, from_side, to_side, reverse]]
                        sub_ref_dict = {to_ref : ls}
                        directed_contact_breakpoint_dict[from_ref][from_pos].update(sub_ref_dict)
                    else:
                        found = False
                        for breakpoint in directed_contact_breakpoint_dict[from_ref][from_pos][to_ref]:
                            if breakpoint[0] == to_pos and breakpoint[-1]  == reverse and breakpoint[1] == from_side and breakpoint[2] == to_side :
                                found = True
                                break
                        if found is False:
                            directed_contact_breakpoint_dict[from_ref][from_pos][to_ref].append([to_pos, from_side, to_side, reverse])
    return directed_contact_breakpoint_dict


def updateBreakpoints(all_seq_dict, candidate_reference_dict, directed_contact_breakpoint_dict, ConnectRefOfEachRefDict):
    bkp_change_map = {}
    junction_edge_list = []
    for seq_name in candidate_reference_dict:
        ref_seq = candidate_reference_dict[seq_name]
        harbor_name = all_seq_dict[seq_name][0][0]
        seq_start_pos = 0
        seq_end_pos = -1
        for seg in ref_seq:
            seg_ref_name = seg[0]
            seg_start_pos = seg[1]
            seg_end_pos = seg[2]
            seq_start_pos = seq_end_pos + 1
            seq_end_pos = seq_start_pos + seg[2] - seg[1]
            if seg[3] == '+':
                seg_reverse = False
            else:
                seg_reverse = True
            if seg_ref_name in directed_contact_breakpoint_dict:
                for from_pos in directed_contact_breakpoint_dict[seg_ref_name]:
                    if from_pos > seg_start_pos and from_pos < seg_end_pos:
                        if seg_reverse is False:#forward insert
                            change_pos = seq_start_pos + from_pos - seg_start_pos
                        else:
                            change_pos = seq_start_pos + (seg_end_pos - from_pos)
                        if seg_ref_name not in bkp_change_map:
                            bkp_change_map.update({seg_ref_name : {from_pos : [[seq_name, change_pos, seg_reverse]]}})
                        else:
                            if from_pos not in bkp_change_map[seg_ref_name]:
                                bkp_change_map[seg_ref_name].update({from_pos : [[seq_name, change_pos, seg_reverse]]})
                            else:
                                if [seq_name, change_pos, seg_reverse] not in bkp_change_map[seg_ref_name][from_pos]:
                                    bkp_change_map[seg_ref_name][from_pos].append([seq_name, change_pos, seg_reverse])                                
            for from_ref in directed_contact_breakpoint_dict:
                for from_pos in directed_contact_breakpoint_dict[from_ref]:
                    if seg_ref_name in directed_contact_breakpoint_dict[from_ref][from_pos]:
                        for to_bkp in directed_contact_breakpoint_dict[from_ref][from_pos][seg_ref_name]:
                            to_pos = to_bkp[0]
                            if to_pos > seg_start_pos and to_pos < seg_end_pos:
                                if seg_reverse is False:
                                    change_pos = seq_start_pos + to_pos - seg_start_pos
                                else:
                                    change_pos = seq_start_pos + (seg_end_pos - to_pos)
                                if seg_ref_name not in bkp_change_map:
                                    bkp_change_map.update({seg_ref_name : {to_pos : [[seq_name, change_pos, seg_reverse]]}})
                                else:
                                    if to_pos not in bkp_change_map[seg_ref_name]:
                                        bkp_change_map[seg_ref_name].update({to_pos : [[seq_name, change_pos, seg_reverse]]})
                                    else:
                                        if [seq_name, change_pos, seg_reverse] not in bkp_change_map[seg_ref_name][to_pos]:
                                            bkp_change_map[seg_ref_name][to_pos].append([seq_name, change_pos, seg_reverse])
    for from_ref in directed_contact_breakpoint_dict:
        for from_pos in directed_contact_breakpoint_dict[from_ref]:
            for to_ref in directed_contact_breakpoint_dict[from_ref][from_pos]:
                for bkp in directed_contact_breakpoint_dict[from_ref][from_pos][to_ref]:
                    to_pos = bkp[0]
                    from_side = bkp[1]
                    to_side = bkp[2]
                    found = False
                    for seq_name in ConnectRefOfEachRefDict:
                        for hgt_event in ConnectRefOfEachRefDict[seq_name]:
                            if hgt_event[0]['from_ref']==from_ref and hgt_event[0]['to_ref']==to_ref and hgt_event[0]['from_pos']==from_pos and hgt_event[0]['to_pos']==to_pos and hgt_event[0]['from_side']==from_side and hgt_event[0]['to_side']==to_side and hgt_event[0]['reverse']==bkp[3]:
                                found = True
                                break
                            if hgt_event[1]['from_ref']==from_ref and hgt_event[1]['to_ref']==to_ref and hgt_event[1]['from_pos']==from_pos and hgt_event[1]['to_pos']==to_pos and hgt_event[1]['from_side']==from_side and hgt_event[1]['to_side']==to_side and hgt_event[1]['reverse']==bkp[3]:
                                found = True
                                break
                            if bkp[3] == True:
                                if hgt_event[0]['to_ref']==from_ref and hgt_event[0]['from_ref']==to_ref and hgt_event[0]['to_pos']==from_pos and hgt_event[0]['from_pos']==to_pos and hgt_event[0]['to_side']==from_side and hgt_event[0]['from_side']==to_side and hgt_event[0]['reverse']==bkp[3]:
                                    found = True
                                    break
                                if hgt_event[1]['to_ref']==from_ref and hgt_event[1]['from_ref']==to_ref and hgt_event[1]['to_pos']==from_pos and hgt_event[1]['from_pos']==to_pos and hgt_event[1]['to_side']==from_side and hgt_event[1]['from_side']==to_side and hgt_event[1]['reverse']==bkp[3]:
                                    found = True
                                    break
                        if found is True:
                            break
                    if found is True:
                        continue
                    if from_ref in bkp_change_map and from_pos in bkp_change_map[from_ref]:
                        for changed_from_bkp_info in bkp_change_map[from_ref][from_pos]:
                            changed_from_ref = changed_from_bkp_info[0]
                            changed_from_pos = changed_from_bkp_info[1]
                            changed_from_seg_reverse = changed_from_bkp_info[2]
                            if changed_from_seg_reverse is True:
                                if from_side == 'right':
                                    changed_from_side = 'left'
                                else:
                                    changed_from_side = 'right'
                            else:
                                changed_from_side = from_side
                            if to_ref in bkp_change_map and to_pos in bkp_change_map[to_ref]:
                                for changed_to_bkp_info in bkp_change_map[to_ref][to_pos]:
                                    changed_to_ref = changed_to_bkp_info[0]
                                    changed_to_pos = changed_to_bkp_info[1]
                                    changed_to_seg_reverse = changed_to_bkp_info[2]
                                    if changed_to_seg_reverse is True:
                                        if to_side == 'right':
                                            changed_to_side = 'left'
                                        else:
                                            changed_to_side = 'right'
                                    else:
                                        changed_to_side = to_side
                                    if changed_from_side == changed_to_side:
                                        changed_reverse = True
                                    else:
                                        changed_reverse = False
                                    junction_edge_list.append([changed_from_ref, changed_from_pos, changed_from_side, changed_to_ref, changed_to_pos, changed_to_side, changed_reverse])
                            else:
                                changed_to_ref = to_ref
                                changed_to_pos = to_pos
                                changed_to_side = to_side
                                if changed_from_side == changed_to_side:
                                    changed_reverse = True
                                else:
                                    changed_reverse = False
                                junction_edge_list.append([changed_from_ref, changed_from_pos, changed_from_side, changed_to_ref, changed_to_pos, changed_to_side, changed_reverse])
                    else:
                        changed_from_ref = from_ref
                        changed_from_pos = from_pos
                        changed_from_side = from_side
                        if to_ref in bkp_change_map and to_pos in bkp_change_map[to_ref]:
                            for changed_to_bkp_info in bkp_change_map[to_ref][to_pos]:
                                changed_to_ref = changed_to_bkp_info[0]
                                changed_to_pos = changed_to_bkp_info[1]
                                changed_to_seg_reverse = changed_to_bkp_info[2]
                                if changed_to_seg_reverse is True:
                                    if to_side == 'right':
                                        changed_to_side = 'left'
                                    else:
                                        changed_to_side = 'right'
                                else:
                                    changed_to_side = to_side
                                if changed_from_side == changed_to_side:
                                    changed_reverse = True
                                else:
                                    changed_reverse = False
                                junction_edge_list.append([changed_from_ref, changed_from_pos, changed_from_side, changed_to_ref, changed_to_pos, changed_to_side, changed_reverse])
                        else:
                            changed_to_ref = to_ref
                            changed_to_pos = to_pos
                            changed_to_side = to_side
                            if changed_from_side == changed_to_side:
                                changed_reverse = True
                            else:
                                changed_reverse = False
                            junction_edge_list.append([changed_from_ref, changed_from_pos, changed_from_side, changed_to_ref, changed_to_pos, changed_to_side, changed_reverse])
    return junction_edge_list, bkp_change_map



#seq为待处理序列
def getOriginalComponentOfSeq(seq, all_seq_dict, genome_sequence_coverage):
    ori_seq = []
    for seg in seq:
        seg_ref_name = seg[0]
        seg_start_pos = seg[1]
        seg_end_pos = seg[2]
        seg_reverse = seg[3]
        if seg_ref_name in genome_sequence_coverage:
            ori_seq.append(seg)
            continue
        ori_seq_start_pos = 0
        ori_seq_end_pos = -1
        tmp_seq_region = []
        for ori_seg in all_seq_dict[seg_ref_name]:
            ori_seg_ref_name = ori_seg[0]
            ori_seg_start_pos = ori_seg[1]
            ori_seg_end_pos = ori_seg[2]
            ori_seq_start_pos = ori_seq_end_pos + 1
            ori_seq_end_pos = ori_seq_start_pos + ori_seg[2] - ori_seg[1]
            ori_seg_reverse = ori_seg[3]
            if ori_seq_start_pos >= seg_start_pos and ori_seq_end_pos <= seg_end_pos:
                tmp_seq_region.append(ori_seg)
                continue
            if ori_seq_start_pos == seg_start_pos and ori_seq_end_pos == seg_end_pos:
                break
            if ori_seq_start_pos <= seg_start_pos and ori_seq_end_pos >= seg_end_pos:
                if ori_seg_reverse is '+':
                    tmp_start_pos = ori_seg_start_pos + seg_start_pos - ori_seq_start_pos
                    tmp_end_pos = ori_seg_end_pos - (ori_seq_end_pos - seg_end_pos)
                else:
                    tmp_start_pos = ori_seg_start_pos + ori_seq_end_pos - seg_end_pos
                    tmp_end_pos = ori_seg_end_pos - (seg_start_pos - ori_seq_start_pos)
                tmp_seq_region.append([ori_seg_ref_name, tmp_start_pos, tmp_end_pos, ori_seg_reverse])
                break
            if ori_seq_start_pos > seg_start_pos and ori_seq_start_pos < seg_end_pos and ori_seq_end_pos > seg_end_pos:
                if ori_seg_reverse is '+':
                    tmp_start_pos = ori_seg_start_pos
                    tmp_end_pos = ori_seg_end_pos - (ori_seq_end_pos - seg_end_pos)
                else:
                    tmp_start_pos = ori_seg_start_pos + ori_seq_end_pos - seg_end_pos
                    tmp_end_pos = ori_seg_end_pos
                tmp_seq_region.append([ori_seg_ref_name, tmp_start_pos, tmp_end_pos, ori_seg_reverse])
                continue
            if ori_seq_start_pos < seg_start_pos and ori_seq_end_pos > seg_start_pos and ori_seq_end_pos < seg_end_pos:
                if ori_seg_reverse is '+':
                    tmp_start_pos = ori_seg_start_pos + seg_start_pos - ori_seq_start_pos
                    tmp_end_pos = ori_seg_end_pos
                else:
                    tmp_start_pos = ori_seg_start_pos
                    tmp_end_pos = ori_seg_end_pos - (seg_start_pos - ori_seq_start_pos)
                tmp_seq_region.append([ori_seg_ref_name, tmp_start_pos, tmp_end_pos, ori_seg_reverse])
                continue
        if seg_reverse is '+':
            for tmp_seg in tmp_seq_region:
                ori_seq.append(tmp_seg)
        else:
            for i in range(len(tmp_seq_region)-1, -1, -1):
                tmp_seg = tmp_seq_region[i]
                if tmp_seg[-1] is '+':
                    tmp_seg[-1] = '-'
                else:
                    tmp_seg[-1] = '+'
                ori_seq.append(tmp_seg)
    return ori_seq

            

def getAllBkp(hgt_seq_dict):
    all_bkp_dict = {}
    for seq_name in hgt_seq_dict:
        for seg in hgt_seq_dict[seq_name]:
            seg_ref_name = seg[0]
            seg_start_pos = seg[1]
            seg_end_pos = seg[2]
            if seg_ref_name not in all_bkp_dict:
                all_bkp_dict.update({seg_ref_name:[seg_start_pos,seg_end_pos]})
            else:
                if seg_start_pos not in all_bkp_dict[seg_ref_name]:
                    all_bkp_dict[seg_ref_name].append(seg_start_pos)
                if seg_end_pos not in all_bkp_dict[seg_ref_name]:
                    all_bkp_dict[seg_ref_name].append(seg_end_pos)
    for bkp_ref_name in all_bkp_dict:
        all_bkp_dict[bkp_ref_name].sort()
    return all_bkp_dict


def getHarborNameOfSeq(hgt_seq_dict):
    harbor2seq_dict = {}
    seq2harbor_dict = {}
    for seq_name in hgt_seq_dict:
        seq2harbor_dict.update({seq_name:hgt_seq_dict[seq_name][0][0]})
    harbor_name_list = []
    for seq_name in seq2harbor_dict:
        if seq2harbor_dict[seq_name] not in harbor_name_list:
            harbor_name_list.append(seq2harbor_dict[seq_name])
    for harbor_name in harbor_name_list:
        ls = []
        for seq_name in seq2harbor_dict:
            if seq2harbor_dict[seq_name] == harbor_name:
                ls.append(seq_name)
        harbor2seq_dict.update({harbor_name:ls})
    return seq2harbor_dict, harbor2seq_dict


def getHarborSeq(harbor2seq_dict, all_bkp_dict):
    harbor_seg_dict = {}
    for bkp_ref_name in all_bkp_dict:
        seg_list = []
        for i in range(len(all_bkp_dict[bkp_ref_name]) - 1):
            seg_start_pos = all_bkp_dict[bkp_ref_name][i]
            seg_end_pos = all_bkp_dict[bkp_ref_name][i + 1]
            if seg_end_pos - seg_start_pos < rlen:
                continue
            seg_list.append([seg_start_pos, seg_end_pos])
        harbor_seg_dict.update({bkp_ref_name:seg_list})
    return harbor_seg_dict




def getAllHgtSeg(hgt_seq_dict, seq2harbor_dict):
    hgt_seg_dict = {}
    for seq_name in hgt_seq_dict:
        harbor_ref_name = seq2harbor_dict[seq_name]
        tmp_ls = []
        for i in range(1, len(hgt_seq_dict[seq_name])):
            seg = hgt_seq_dict[seq_name][i]
            if seg[0] != harbor_ref_name:
                tmp_ls.append(i)
        if len(tmp_ls)>0:
            hgt_seg_dict.update({seq_name : tmp_ls})
    return hgt_seg_dict




def mapSeg2Seq(all_bkp_dict, hgt_seq_dict, hgt_seg_dict, seq2harbor_dict):
    seq_cov_dict = {}
    for bkp_ref_name in all_bkp_dict:
        non_hgt_seg_list = []
        for i in range(len(all_bkp_dict[bkp_ref_name]) - 1):
            seg_start_pos = all_bkp_dict[bkp_ref_name][i]
            seg_end_pos = all_bkp_dict[bkp_ref_name][i + 1]
            if seg_end_pos - seg_start_pos < rlen:
                continue
            found = False
            for seq_name in hgt_seg_dict:
                if seq2harbor_dict[seq_name]==bkp_ref_name:
                    continue
                for index in hgt_seg_dict[seq_name]:
                    seg = hgt_seq_dict[seq_name][index]
                    if seg[0] == bkp_ref_name:
                        if seg_start_pos >= seg[1] and seg_end_pos <= seg[2]:
                            found = True
                            break
                if found is True:
                    break
            if found is False:
                non_hgt_seg_list.append([seg_start_pos, seg_end_pos])
        avg_cov = 0
        seg_num = 0
        if len(non_hgt_seg_list) > 0:
            for seg in non_hgt_seg_list:
                tmp_cov = getSegCoverage(bkp_ref_name, seg[0], seg[1])
                if math.isnan(tmp_cov):
                    continue
                avg_cov += tmp_cov
                seg_num += 1
            if seg_num != 0:
                avg_cov = avg_cov/seg_num
        else:
            avg_cov = getRefCoverage(bkp_ref_name)
        if math.isnan(avg_cov) or seg_num == 0:
            avg_cov = getRefCoverage(bkp_ref_name)
        seq_cov_dict.update({bkp_ref_name : avg_cov})
    seg_cov_dict = {}
    seg2seq_dict = {}
    for bkp_ref_name in all_bkp_dict:
        tmp_cov_dict = {}
        tmp_dict = {}
        for i in range(len(all_bkp_dict[bkp_ref_name]) - 1):
            seg_start_pos = all_bkp_dict[bkp_ref_name][i]
            seg_end_pos = all_bkp_dict[bkp_ref_name][i + 1]
            if seg_end_pos - seg_start_pos < rlen:
                continue
            ls = []
            for seq_name in hgt_seg_dict:
                if seq2harbor_dict[seq_name]==bkp_ref_name:
                    continue
                for index in hgt_seg_dict[seq_name]:
                    seg = hgt_seq_dict[seq_name][index]
                    if seg[0] == bkp_ref_name:
                        if seg_start_pos >= seg[1] and seg_end_pos <= seg[2]:
                            ls.append([seq_name, index])
            if len(ls) > 0:
                #avg_cov = getSegCoverage(bkp_ref_name, seg_start_pos, seg_end_pos)
                seg_cov = 0
                for seg_info in ls:
                    seg_cov += seq_cov_dict[seq2harbor_dict[seg_info[0]]]
                #ref_cov = getRefCoverage(bkp_ref_name)
                #if avg_cov < ref_cov:
                #    avg_cov += ref_cov
                found = False
                for index_key in tmp_dict:
                    if tmp_dict[index_key] == ls:
                        tmp_cov_dict[index_key].append(seg_cov)
                        found = True
                        break
                if found is True:
                    continue
                tmp_cov_dict.update({i:[seg_cov]})
                tmp_dict.update({i:ls})
        if len(tmp_dict) == 0:
            continue
        seg2seq_dict.update({bkp_ref_name:tmp_dict})
        seg_cov_dict.update({bkp_ref_name:tmp_cov_dict})
    return seg2seq_dict, seg_cov_dict, seq_cov_dict



def compareList(ls1, ls2):
    if len(ls1) == len(ls2):
        for seg in ls1:
            if seg not in ls2:
                return False
        return True
    return False 


def getJunctionEdge(harbor2seq_dict, seq2harbor_dict, harbor_seg_dict, hgt_seq_dict):
    out_junction_edge_dict = {}
    for harbor_name in harbor2seq_dict:
        tmp_dict = {}
        for i in range(len(harbor_seg_dict[harbor_name])-1):
            seg = harbor_seg_dict[harbor_name][i]
            ls = []
            for seq_name in harbor2seq_dict[harbor_name]:
                for j in range(len(hgt_seq_dict[seq_name])-1):
                    if hgt_seq_dict[seq_name][j][0] == harbor_name and hgt_seq_dict[seq_name][j+1][0] != harbor_name and abs(hgt_seq_dict[seq_name][j][2] - seg[1])<rlen:
                        if [seq_name, j+1] not in ls:
                            ls.append([seq_name, j+1])
            if len(ls)>0:
                tmp_dict.update({i:ls})
        if any(tmp_dict) is True:
            out_junction_edge_dict.update({harbor_name:tmp_dict})
    back_junction_edge_dict = {}
    for harbor_name in harbor2seq_dict:
        tmp_dict = {}
        for i in range(1, len(harbor_seg_dict[harbor_name])):
            seg = harbor_seg_dict[harbor_name][i]
            ls = []
            for seq_name in harbor2seq_dict[harbor_name]:
                for j in range(1, len(hgt_seq_dict[seq_name])):
                    if hgt_seq_dict[seq_name][j][0] == harbor_name and hgt_seq_dict[seq_name][j-1][0] != harbor_name and abs(hgt_seq_dict[seq_name][j][1] - seg[0])<rlen:
                        if [seq_name, j-1] not in ls:
                            ls.append([seq_name, j-1])
            if len(ls)>0:
                found = False
                if harbor_name in out_junction_edge_dict:
                    for index in out_junction_edge_dict[harbor_name]:
                        if compareList(out_junction_edge_dict[harbor_name][index], ls):
                            found = True
                            break
                if found is True:
                    continue
                tmp_dict.update({i:ls})
        if any(tmp_dict) is True:
            back_junction_edge_dict.update({harbor_name:tmp_dict})
    near_junction_edge_dict = {}
    for seq_name in hgt_seq_dict:
        tmp_dict = {}
        for i in range(len(hgt_seq_dict[seq_name])-1):
            if hgt_seq_dict[seq_name][i][0]!=seq2harbor_dict[seq_name] and hgt_seq_dict[seq_name][i+1][0]!=seq2harbor_dict[seq_name] and hgt_seq_dict[seq_name][i][0]!=hgt_seq_dict[seq_name][i+1][0]:
                tmp_dict.update({i:i+1})
        if any(tmp_dict) is True:
            near_junction_edge_dict.update({seq_name:tmp_dict})
    return out_junction_edge_dict, back_junction_edge_dict, near_junction_edge_dict


def getJunctionGraph(harbor2seq_dict, harbor_seg_dict, hgt_seq_dict):
    junction_graph_dict = {}
    for harbor_name in harbor2seq_dict:
        tmp_dict = {}
        for seq_name in harbor2seq_dict[harbor_name]:
            split_seg_list = []
            for i in range(len(hgt_seq_dict[seq_name])):
                if hgt_seq_dict[seq_name][i][0] == harbor_name:
                    for j in range(0, len(harbor_seg_dict[harbor_name])):
                        seg = harbor_seg_dict[harbor_name][j]
                        if seg[0] > hgt_seq_dict[seq_name][i][1] - rlen and seg[1] < hgt_seq_dict[seq_name][i][2] + rlen:
                            split_seg_list.append([harbor_name, j])
                else:
                    split_seg_list.append([seq_name, i])
            for i in range(len(split_seg_list)-1):
                if split_seg_list[i][0] not in tmp_dict:
                    tmp_dict.update({split_seg_list[i][0]:{split_seg_list[i][1]:{split_seg_list[i+1][0]:[split_seg_list[i+1][1]]}}})
                else:
                    if split_seg_list[i][1] not in tmp_dict[split_seg_list[i][0]]:
                        tmp_dict[split_seg_list[i][0]].update({split_seg_list[i][1]:{split_seg_list[i+1][0]:[split_seg_list[i+1][1]]}})
                    else:
                        if split_seg_list[i+1][0] not in tmp_dict[split_seg_list[i][0]][split_seg_list[i][1]]:
                            tmp_dict[split_seg_list[i][0]][split_seg_list[i][1]].update({split_seg_list[i+1][0]:[split_seg_list[i+1][1]]})
                        else:
                            if split_seg_list[i+1][1] not in tmp_dict[split_seg_list[i][0]][split_seg_list[i][1]][split_seg_list[i+1][0]]:
                                tmp_dict[split_seg_list[i][0]][split_seg_list[i][1]][split_seg_list[i+1][0]].append(split_seg_list[i+1][1])
        if any(tmp_dict) is True:
            junction_graph_dict.update({harbor_name:tmp_dict})
    return junction_graph_dict




def getAllPossibleSeq(junction_graph_dict, junction_seg_cov_dict, harbor_seg_dict, hgt_seq_dict):
    all_possible_hgt_dict = {}
    all_possible_hgt_freq_dict = {}
    for harbor_name in junction_graph_dict:
        hgt_freq_list = []
        hgt_seq_list = []
        count = 0
        while junction_seg_cov_dict[harbor_name][harbor_name][0] > 0:
            hgt_seq = []
            junction_seg_cov_dict[harbor_name][harbor_name][0] -= 1
            tmp_seg = [harbor_name, harbor_seg_dict[harbor_name][0][0], harbor_seg_dict[harbor_name][0][1], '+']
            hgt_seq.append(tmp_seg)
            seg_seq_name = harbor_name
            seg_index = 0
            found = True
            while seg_seq_name != harbor_name or seg_index != len(harbor_seg_dict[harbor_name]) - 1:
                candidate_seg_list = []
                for next_seg_seq_name in junction_graph_dict[harbor_name][seg_seq_name][seg_index]:
                    for next_seg_index in junction_graph_dict[harbor_name][seg_seq_name][seg_index][next_seg_seq_name]:
                        if junction_seg_cov_dict[harbor_name][next_seg_seq_name][next_seg_index] >0:
                            candidate_seg_list.append([next_seg_seq_name, next_seg_index])
                if len(candidate_seg_list) == 0:
                    found = False
                    break
                if len(candidate_seg_list) != 1:
                    random.shuffle(candidate_seg_list)
                seg_seq_name = candidate_seg_list[0][0]
                seg_index = candidate_seg_list[0][1]
                if seg_seq_name!=harbor_name:
                    hgt_seq.append(hgt_seq_dict[seg_seq_name][seg_index])
                else:
                    tmp_seg = [harbor_name, harbor_seg_dict[harbor_name][seg_index][0], harbor_seg_dict[harbor_name][seg_index][1], '+']
                    hgt_seq.append(tmp_seg)
                junction_seg_cov_dict[harbor_name][seg_seq_name][seg_index] -= 1
            if found is False:
                break
            if hgt_seq not in hgt_seq_list:
                hgt_seq_list.append(hgt_seq)
                hgt_freq_list.append(1)
            else:
                seq_index = hgt_seq_list.index(hgt_seq)
                hgt_freq_list[seq_index] += 1
            count += 1
        for i in range(len(hgt_freq_list)):
            hgt_freq_list[i] = hgt_freq_list[i]/count
        all_possible_hgt_dict.update({harbor_name:hgt_seq_list})
        all_possible_hgt_freq_dict.update({harbor_name:hgt_freq_list})
    return all_possible_hgt_dict, all_possible_hgt_freq_dict
        


def getJunctionSegCov(seg_var_dict, seq_cov_dict, harbor2seq_dict, harbor_seg_dict):
    junction_seg_cov_dict = {}
    for harbor_name in harbor2seq_dict:
        tmp_dict = {}
        for j in range(len(harbor_seg_dict[harbor_name])):
            if harbor_name not in tmp_dict:
                tmp_dict.update({harbor_name:{j:seq_cov_dict[harbor_name]}})
            else:
                tmp_dict[harbor_name].update({j:seq_cov_dict[harbor_name]})
        for seq_name in harbor2seq_dict[harbor_name]:
            for j in seg_var_dict[seq_name]:
                if seq_name not in tmp_dict:
                    tmp_dict.update({seq_name:{j:seg_var_dict[seq_name][j].x}})
                else:
                    tmp_dict[seq_name].update({j:seg_var_dict[seq_name][j].x})
        if any(tmp_dict) is True:
            junction_seg_cov_dict.update({harbor_name:tmp_dict})
    return junction_seg_cov_dict



def getCoverageOfSeq(seg2seq_dict, harbor2seq_dict, harbor_seg_dict, seg_cov_dict, seq_cov_dict, hgt_seg_dict, out_junction_edge_dict, back_junction_edge_dict, near_junction_edge_dict):
    m = Model()
    seg_var_dict = {}
    for seq_name in hgt_seg_dict:
        tmp_dict = {}
        for index in hgt_seg_dict[seq_name]:
            seg_var = m.addVar(lb=0, vtype=GRB.CONTINUOUS, name = seq_name+'_'+str(index))
            tmp_dict.update({index:seg_var})
        seg_var_dict.update({seq_name:tmp_dict})
    slack_list_1 = []
    for bkp_ref_name in seg2seq_dict:
        for seg_index in seg2seq_dict[bkp_ref_name]:
            seg_var_list = []
            for seg_info in seg2seq_dict[bkp_ref_name][seg_index]:
                seg_var_list.append(seg_var_dict[seg_info[0]][seg_info[1]])
            slack = m.addVar(lb=0, vtype=GRB.CONTINUOUS)
            seg_coverage = 0
            for cov in seg_cov_dict[bkp_ref_name][seg_index]:
                seg_coverage += cov
            seg_coverage = seg_coverage/len(seg_cov_dict[bkp_ref_name][seg_index])
            m.addConstr(quicksum(seg_var_list) <= seg_coverage + slack)
            m.addConstr(quicksum(seg_var_list) >= seg_coverage - slack)
            slack_list_1.append(slack)
    for ref_key in out_junction_edge_dict:
        for seg_index in out_junction_edge_dict[ref_key]:
            ref_cov = seq_cov_dict[ref_key]
            seg_var_list = []
            for seg_info in out_junction_edge_dict[ref_key][seg_index]:
                seg_var_list.append(seg_var_dict[seg_info[0]][seg_info[1]])
            m.addConstr(quicksum(seg_var_list) == ref_cov)
    for ref_key in back_junction_edge_dict:
        for seg_index in back_junction_edge_dict[ref_key]:
            ref_cov = seq_cov_dict[ref_key]
            seg_var_list = []
            for seg_info in back_junction_edge_dict[ref_key][seg_index]:
                seg_var_list.append(seg_var_dict[seg_info[0]][seg_info[1]])
            m.addConstr(quicksum(seg_var_list) == ref_cov)
    for seq_name in near_junction_edge_dict:
        for index in near_junction_edge_dict[seq_name]:
            m.addConstr(seg_var_dict[seq_name][index] == seg_var_dict[seq_name][near_junction_edge_dict[seq_name][index]])
    m.setObjective(quicksum(slack_list_1), GRB.MINIMIZE)
    m.update()
    m.optimize()
    junction_seg_cov_dict = getJunctionSegCov(seg_var_dict, seq_cov_dict, harbor2seq_dict, harbor_seg_dict)
    return junction_seg_cov_dict



def get_bkp_dict(directed_contact_breakpoint_dict):
    from_ref_list = []
    to_ref_list = []
    ref_bkp_dict = {}
    for from_ref in directed_contact_breakpoint_dict:
        if from_ref not in from_ref_list:
            from_ref_list.append(from_ref)
        if from_ref not in ref_bkp_dict:
            tmp_dict = {from_ref : list(directed_contact_breakpoint_dict[from_ref].keys())}
            ref_bkp_dict.update(tmp_dict)
    for from_ref in directed_contact_breakpoint_dict:
        for from_pos in directed_contact_breakpoint_dict[from_ref]:
            for to_ref in directed_contact_breakpoint_dict[from_ref][from_pos]:
                if to_ref not in to_ref_list:
                    to_ref_list.append(to_ref)
                for breakpoint_info in directed_contact_breakpoint_dict[from_ref][from_pos][to_ref]:
                    if to_ref not in ref_bkp_dict:
                        tmp_dict = {to_ref : [breakpoint_info[0]]}
                        ref_bkp_dict.update(tmp_dict)
                    else:
                        if breakpoint_info[0] not in ref_bkp_dict[to_ref]:
                            ref_bkp_dict[to_ref].append(breakpoint_info[0])
    single_bkp_ref_list = []
    for ref in ref_bkp_dict:
        if len(ref_bkp_dict[ref]) == 1:
            single_bkp_ref_list.append(ref)
    for ref in single_bkp_ref_list:
        del ref_bkp_dict[ref]
    single_direction_ref_list = []
    for ref in ref_bkp_dict:
        if ref not in from_ref_list and ref in to_ref_list:
            single_direction_ref_list.append(ref)
        if ref in from_ref_list and ref not in to_ref_list:
            single_direction_ref_list.append(ref)
    for ref in single_direction_ref_list:
        del ref_bkp_dict[ref]
    for ref_name in ref_bkp_dict:
        candidate_bkp_list = ref_bkp_dict[ref_name]
        if ref_name in directed_contact_breakpoint_dict:
            if abs(0 - min(ref_bkp_dict[ref_name])) > 100:
                ref_bkp_dict[ref_name].append(0)
        ref_bkp_dict[ref_name].sort()
    return ref_bkp_dict




def getDiffRef(currentVisitedRefList):
    visited_ref = []
    for junction_edge in currentVisitedRefList:
        if junction_edge['from_ref'] not in visited_ref:
            visited_ref.append(junction_edge['from_ref'])
        if junction_edge['to_ref'] not in visited_ref:
            visited_ref.append(junction_edge['to_ref'])
    return len(visited_ref)



def getOriSegInfo(seq_name, pos):
    seq_start_pos = 0
    seq_end_pos = -1
    if seq_name in all_seq_dict:
        for seg in all_seq_dict[seq_name]:
            seg_ref_name = seg[0]
            seq_start_pos = seq_end_pos + 1
            seq_end_pos = seq_start_pos + seg[2] - seg[1]
            if pos >= seq_start_pos + rlen and pos <= seq_end_pos + rlen:
                return seg_ref_name
    else:
        return seq_name


def BFS_modify(from_ref, directed_contact_breakpoint_dict, ref_num_threshold, currentVisitedRefList, target_ref, connectedRefList, ref_bkp_dict):
    num_diff_ref = getDiffRef(currentVisitedRefList)
    if len(currentVisitedRefList) > 1:
        return connectedRefList
    if num_diff_ref <= ref_num_threshold:
        for from_pos in directed_contact_breakpoint_dict[from_ref]:
            for to_ref in directed_contact_breakpoint_dict[from_ref][from_pos]:
                if from_ref == to_ref:
                    continue
                if len(currentVisitedRefList) == 1:
                    if currentVisitedRefList[0]['from_ref'] != to_ref:
                        continue
                if to_ref not in ref_bkp_dict:
                    continue
                for bkp in directed_contact_breakpoint_dict[from_ref][from_pos][to_ref]:
                    to_pos = bkp[0]
                    from_side = bkp[1]
                    to_side = bkp[2]
                    reverse = bkp[3]
                    if len(currentVisitedRefList) == 1:
                        if to_ref != currentVisitedRefList[0]['from_ref']:
                            continue
                        if to_ref in all_seq_dict:
                            harbor_name = all_seq_dict[to_ref][0][0]
                        else:
                            harbor_name = to_ref
                        from_seg_name = getOriSegInfo(from_ref, from_pos)
                        to_seg_name = getOriSegInfo(to_ref, to_pos)
                        if to_seg_name!=harbor_name and from_seg_name!=harbor_name:
                            continue
                        previous_to_side = currentVisitedRefList[0]['to_side']
                        previous_reverse = currentVisitedRefList[0]['reverse']
                        if previous_to_side == 'right':
                            if from_side != 'left':
                                continue
                        else:
                            if from_side != 'right':
                                continue
                        if reverse != previous_reverse:
                            continue
                        if to_pos < currentVisitedRefList[0]['from_pos']:
                            continue
                        insert_length = abs(currentVisitedRefList[0]['to_pos'] - from_pos)
                        bkp_length = abs(currentVisitedRefList[0]['from_pos'] - to_pos)
                        if bkp_length > 1000 or insert_length < rlen:
                            continue
                        if reverse is False and currentVisitedRefList[0]['to_pos'] > from_pos:
                            continue
                        harbor_ref_name = to_ref
                        harbor_ref_length = sim_seq_length_dict[harbor_ref_name]
                        #if insert_length / harbor_ref_length > 0.4:
                        #    continue
                        #if harbor_ref_length + insert_length > 20000000:
                        #    continue
                    if to_ref == target_ref:
                        currentVisitedRefList.append({'from_ref':from_ref, 'to_ref':to_ref, 'from_pos':from_pos, 'to_pos':to_pos, 'from_side':from_side, 'to_side':to_side, 'reverse':reverse})#from_ref, to_ref, from_pos, tp_pos
                        tmp_ls = copy.deepcopy(currentVisitedRefList)
                        connectedRefList.append(tmp_ls)
                        currentVisitedRefList.pop()
                        return connectedRefList
                    else:
                        currentVisitedRefList.append({'from_ref':from_ref, 'to_ref':to_ref, 'from_pos':from_pos, 'to_pos':to_pos, 'from_side':from_side, 'to_side':to_side, 'reverse':reverse})
                        connectedRefList = BFS_modify(to_ref, directed_contact_breakpoint_dict, ref_num_threshold, currentVisitedRefList, target_ref, connectedRefList, ref_bkp_dict)
                        currentVisitedRefList.pop()
    return connectedRefList


def getConnectRefOfEachRef(ref_bkp_dict, directed_contact_breakpoint_dict):
    ConnectRefOfEachRefDict = {}
    ref_num_threshold = 2
    for from_ref in ref_bkp_dict:
        if ref_bkp_dict[from_ref][0] != 0:
            continue
        if from_ref in directed_contact_breakpoint_dict:
            connectedRefList = []
            currentVisitedRefList = []
            connectedRefList = BFS_modify(from_ref, directed_contact_breakpoint_dict, ref_num_threshold, currentVisitedRefList, from_ref, connectedRefList, ref_bkp_dict)
            if len(connectedRefList) == 0:
                continue
            tmp_dict = {from_ref : connectedRefList}
            ConnectRefOfEachRefDict.update(tmp_dict)
    return ConnectRefOfEachRefDict



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


def filter_windows(X, Y, window = 1000, inc_threshold = 0.60, mdiff = float(8), wdiff = float(1.5)):
    """
    filter low and high coverage bins based on running average
    and difference from median
    """
    filtered = [[], []] # store x and y
    weights = np.ones(window)/window
    med = np.median(Y)
    avgs = signal.fftconvolve(Y, weights, 'same').tolist()
    for xi, avi, yi in zip(X, avgs, Y):
        # skip if zero
        if yi <= 0 or avi <= 0 or med <= 0:
            continue
        # skip if >wdiff different from average including 1000 coverage windows
        if abs(float(max([yi, avi]))/float(min([yi, avi]))) > wdiff:
            continue
        # skip if >mdiff different from median
        if abs(float(max([yi, med]))/float(min([yi, med]))) > mdiff:
            continue
        filtered[0].append(xi)
        filtered[1].append(yi)
    # > inc_threashold windows must remain
    if float(len(filtered[0])) / float(len(X)) < inc_threshold:
        return False
    return filtered



def coverage_windows(cov, window = 10000, slide = 100):
    """
    sliding window smoothing of coverage data
    """
    # calculate coverage windows for sample
    weights = np.ones(window)
    if len(cov) < len(weights):
        filtered_cov = False
        return filtered_cov
    windows = [[], []] # x and y values for windows
    i = 0
    for c in signal.fftconvolve(cov, weights, 'valid').tolist()[0::slide]:
        windows[0].append(i)
        windows[1].append(c/window)
        i += slide
    # filter high and low coverage windows
    filtered = filter_windows(windows[0], windows[1])
    if filtered is False:
        filtered_cov = False
    else:
        # log transform
        filtered_cov = [filtered[0], filtered[1]]
    return filtered_cov


def median_filter(y, window = 999):
    """
    return median filtered data
    """
    return ndimage.filters.median_filter(y, size = window, mode = 'reflect')


def get_median_filtered_coverage(region_cov):
    filtered_cov = coverage_windows(region_cov, window = 10000, slide = 100)
    if filtered_cov is False:
        return median_filter(region_cov)
    else:
        return median_filter(filtered_cov[1])



def getTargetSeg(ref_name, out_back_info, candidate_reference, start_pos):
    interact_ref_name = out_back_info['interact_ref_name']
    current_out_pos = out_back_info['out_pos']
    current_back_pos = out_back_info['back_pos']
    target_pos = out_back_info['target_pos']
    back_from_pos = out_back_info['back_from_pos']
    #if abs(back_from_pos - target_pos) > 10000:
    if start_pos < min(current_out_pos, current_back_pos):
        candidate_reference.append([ref_name, start_pos, min(current_out_pos, current_back_pos), '+'])
    current_max_pos = max(current_out_pos, current_back_pos)
    if out_back_info['out_side'] == 'right':#from ref is '+'
        if out_back_info['out_reverse'] == False:
            if back_from_pos > target_pos:
                if [interact_ref_name, target_pos, back_from_pos, '+'] != candidate_reference[-1]:
                    candidate_reference.append([interact_ref_name, target_pos, back_from_pos, '+'])
        else:
            if back_from_pos < target_pos:
                if [interact_ref_name, back_from_pos, target_pos, '-'] != candidate_reference[-1]:
                    candidate_reference.append([interact_ref_name, back_from_pos, target_pos, '-'])
    else:#from_ref is '-'
        if out_back_info['out_reverse'] == True:
            if back_from_pos > target_pos:
                if [interact_ref_name, target_pos, back_from_pos, '-'] != candidate_reference[-1]:
                    candidate_reference.append([interact_ref_name, target_pos, back_from_pos, '-'])
    if current_out_pos == current_back_pos:
        current_max_pos += 1
    return candidate_reference, current_max_pos


def getRefNum(candidate_reference):
    ref_name_list = []
    for segment in candidate_reference:
        if segment[0] not in ref_name_list:
            ref_name_list.append(segment[0])
    return len(ref_name_list)


def getSeqLength(seq):
    length = 0
    for seg in seq:
        length += seg[2] - seg[1] + 1
    return length


def getAssembleReference(ConnectRefOfEachRefDict, order):
    candidate_reference_dict = {}
    count = 0
    for ref_name in ConnectRefOfEachRefDict:
        edge_path_list = ConnectRefOfEachRefDict[ref_name]
        out_back_info_list = []
        for edge_path in edge_path_list:
            if edge_path[0]['from_ref'] == ref_name and edge_path[1]['to_ref'] == ref_name:
                out_pos = edge_path[0]['from_pos']
                out_side = edge_path[0]['from_side']
                out_reverse = edge_path[0]['reverse']
                interact_ref_name = edge_path[0]['to_ref']
                target_pos = edge_path[0]['to_pos']
                back_pos = edge_path[1]['to_pos']
                back_side = edge_path[1]['to_side']
                back_reverse = edge_path[1]['reverse']
                back_from_pos = edge_path[1]['from_pos']
                if back_from_pos > sim_seq_length_dict[interact_ref_name]:
                    back_from_pos = sim_seq_length_dict[interact_ref_name]
                if target_pos > sim_seq_length_dict[interact_ref_name]:
                    target_pos = sim_seq_length_dict[interact_ref_name]
                if out_pos > sim_seq_length_dict[ref_name]:
                    out_pos = sim_seq_length_dict[ref_name]
                if back_pos > sim_seq_length_dict[ref_name]:
                    back_pos = sim_seq_length_dict[ref_name]
                out_back_info = {'out_pos' : out_pos, 'out_side' : out_side, 'out_reverse' : out_reverse, 'interact_ref_name' : interact_ref_name, 'target_pos' : target_pos, 'back_pos' : back_pos, 'back_side' : back_side, 'back_reverse' : back_reverse, 'back_from_pos' : back_from_pos, 'visited' : False}
                #out_back_info = [out_pos, out_side, out_reverse, target_pos, back_pos, back_side, back_reverse, back_from_pos, interact_ref_name, False]
                if out_back_info not in out_back_info_list:
                    out_back_info_list.append(out_back_info)
        out_back_info_list.sort(key=lambda ele:ele['out_pos'])
        #out_back_info_list.sort()
        candidate_reference_list = []
        for i in range(0, len(out_back_info_list)):
            out_back_info = out_back_info_list[i]
            if out_back_info['visited'] == False:
                candidate_reference = []
                if out_back_info['out_side'] == 'left' and out_back_info['out_pos'] < out_back_info['back_pos']:
                    continue
                if out_back_info['out_side'] == 'right' and out_back_info['out_pos'] > out_back_info['back_pos']:
                    continue
                candidate_reference, current_max_pos = getTargetSeg(ref_name, out_back_info, candidate_reference, 0)
                out_back_info['visited'] = True
                if i != len(out_back_info_list) - 1:
                    for j in range(i + 1, len(out_back_info_list)):
                        if out_back_info_list[j]['out_pos'] < current_max_pos + 100 or out_back_info_list[j]['back_pos'] < current_max_pos + 100:
                            continue
                        out_back_info = out_back_info_list[j]
                        if out_back_info['visited'] == True:
                            continue
                        if out_back_info['out_side'] == 'left' and out_back_info['out_pos'] < out_back_info['back_pos']:
                            continue
                        if out_back_info['out_side'] == 'right' and out_back_info['out_pos'] > out_back_info['back_pos']:
                            continue
                        candidate_reference, current_max_pos = getTargetSeg(ref_name, out_back_info, candidate_reference, current_max_pos)
                        out_back_info['visited'] = True
                        if getSeqLength(candidate_reference) > 20000000:
                            break
                if len(candidate_reference) != 0:
                    if current_max_pos < sim_seq_length_dict[ref_name]:
                        candidate_reference.append([ref_name, current_max_pos, sim_seq_length_dict[ref_name], '+'])
                        if candidate_reference not in candidate_reference_list and getSeqLength(candidate_reference) < 20000000:
                            RefNum = getRefNum(candidate_reference)
                            if RefNum > 1:
                                hgt_seq_name = str(order) + '_' + str(count)
                                tmp_dict = {hgt_seq_name : candidate_reference}
                                candidate_reference_dict.update(tmp_dict)
                                count += 1
    return candidate_reference_dict


def calRefCoverage(ref_name):
    ref_cov_array = np.asarray(list(genome_sequence_coverage[ref_name].values()))
    filtered_coverage = get_median_filtered_coverage(ref_cov_array)
    average_coverage = np.average(filtered_coverage)
    return average_coverage


def getRefCoverage(ref_name):
    if ref_name in ref_coverage_dict:
        ref_coverage = ref_coverage_dict[ref_name]
    else:
        ref_coverage = calRefCoverage(ref_name)
        tmp_dict = {ref_name : ref_coverage}
        ref_coverage_dict.update(tmp_dict)
    return ref_coverage


def consumeAndGetRefCoverage(ref_name, CostRefCoverage):
    if ref_name not in remained_coverage_dict:
        ref_coverage = getRefCoverage(ref_name)
    else:
        ref_coverage = remained_coverage_dict[ref_name]
    remained_coverage_dict.update({ref_name : max(ref_coverage - CostRefCoverage, 0)})
    return ref_coverage



def calSegCoverage(ref_name, left_pos, right_pos):
    tmp_coverage_list = []
    for position in range(left_pos, right_pos):
        if position in genome_sequence_coverage[ref_name]:
            tmp_coverage_list.append(genome_sequence_coverage[ref_name][position])
    tmp_coverage_array = np.asarray(tmp_coverage_list)
    filtered_coverage = get_median_filtered_coverage(tmp_coverage_array)
    average_coverage = np.average(filtered_coverage)
    return average_coverage


def getSegCoverage(ref_name, left_pos, right_pos):
    if ref_name in seg_coverage_dict and left_pos in seg_coverage_dict[ref_name] and right_pos in seg_coverage_dict[ref_name][left_pos]:
        return seg_coverage_dict[ref_name][left_pos][right_pos]
    else:
        avg_cov = calSegCoverage(ref_name, left_pos, right_pos)
        if ref_name not in seg_coverage_dict:
            tmp_dict = {ref_name : {left_pos : {right_pos : avg_cov}}}
            seg_coverage_dict.update(tmp_dict)
        else:
            if left_pos not in seg_coverage_dict[ref_name]:
                tmp_dict = {left_pos : {right_pos : avg_cov}}
                seg_coverage_dict[ref_name].update(tmp_dict)
            else:
                if right_pos not in seg_coverage_dict[ref_name][left_pos]:
                    tmp_dict = {right_pos : avg_cov}
                    seg_coverage_dict[ref_name][left_pos].update(tmp_dict)
        return avg_cov



def readTruthSeq(file_name):
    true_hgt_dict = {}
    true_hgt_seg_dict = {}
    tmp_hgt_dict = {}
    fi = open(file_name, "r")
    seq_name = ''
    ls=[]
    while 1:
        buf = fi.readline()
        buf = buf.strip('\n')
        if not buf:
            break
        buf = re.split(r'[:;,\s]\s*', buf)
        if len(buf) == 1:
            if len(ls) > 0:
                tmp_hgt_dict.update({seq_name : ls})
                ls = []
            seq_name = buf[0]
            continue
        ls.append([buf[0], int(buf[1]), int(buf[2]), buf[3]])
    for ref_key in tmp_hgt_dict:
        ls = tmp_hgt_dict[ref_key]
        receptor_name = ls[0][0]
        hgt_seg_list = []
        i = 0
        while i < len(ls)-1:
            if ls[i+1][0] != receptor_name:
                hgt_seg = []
                for j in range(i, len(ls)):
                    if ls[j][0] == receptor_name and ls[j-1][0] != receptor_name and len(hgt_seg) > 0:
                        hgt_seg.append(ls[j])
                        hgt_seg_list.append(hgt_seg)
                        i = j
                        break
                    hgt_seg.append(ls[j])
            else:
                i = i + 1
        if len(hgt_seg_list) > 0:
            true_hgt_seg_dict.update({receptor_name : hgt_seg_list})
    tmp_hgt_dict.update({seq_name : ls})
    for seq_name in tmp_hgt_dict:
        harbor_name = tmp_hgt_dict[seq_name][0][0]
        true_hgt_dict.update({harbor_name:tmp_hgt_dict[seq_name]})
    return true_hgt_dict, true_hgt_seg_dict



def split_segment(trueHGT, simHgtList):
    delta = 20
    change_map_dict = {}
    breakpoint_dict = {}
    seq_list = []
    for seq in simHgtList:
        seq_list.append(seq)
    seq_list.append(trueHGT)
    for seq in seq_list:
        for seg in seq:
            if seg[0] not in change_map_dict:
                map_dict = {seg[1] : seg[1]}
                tmp_dict = {seg[0] : map_dict}
                change_map_dict.update(tmp_dict)
                if seg[2] not in change_map_dict[seg[0]]:
                    found = False
                    for pos_key in change_map_dict[seg[0]]:
                        if abs(change_map_dict[seg[0]][pos_key] - seg[2]) < delta:
                            map_dict = {seg[2] : change_map_dict[seg[0]][pos_key]}
                            change_map_dict[seg[0]].update(map_dict)
                            found = True
                            break
                    if found is False:
                        map_dict = {seg[2] : seg[2]}
                        change_map_dict[seg[0]].update(map_dict)
            else:
                if seg[1] not in change_map_dict[seg[0]]:
                    found = False
                    for pos_key in change_map_dict[seg[0]]:
                        if abs(change_map_dict[seg[0]][pos_key] - seg[1]) < delta:
                            map_dict = {seg[1] : change_map_dict[seg[0]][pos_key]}
                            change_map_dict[seg[0]].update(map_dict)
                            found = True
                            break
                    if found is False:
                        map_dict = {seg[1] : seg[1]}
                        change_map_dict[seg[0]].update(map_dict)
                if seg[2] not in change_map_dict[seg[0]]:
                    found = False
                    for pos_key in change_map_dict[seg[0]]:
                        if abs(change_map_dict[seg[0]][pos_key] - seg[2]) < delta:
                            map_dict = {seg[2] : change_map_dict[seg[0]][pos_key]}
                            change_map_dict[seg[0]].update(map_dict)
                            found = True
                            break
                    if found is False:
                        map_dict = {seg[2] : seg[2]}
                        change_map_dict[seg[0]].update(map_dict)
    for seq in seq_list:
        for seg in seq:
            seg[1] = change_map_dict[seg[0]][seg[1]]
            seg[2] = change_map_dict[seg[0]][seg[2]]
    for seq in seq_list:
        for seg in seq:
            if seg[0] not in breakpoint_dict:
                tmp_dict = {seg[0] : [seg[1], seg[2]]}
                breakpoint_dict.update(tmp_dict)
            else:
                found = False
                for bkp in breakpoint_dict[seg[0]]:
                    if abs(bkp - seg[1]) < delta:
                        found = True
                        break
                if found is False:
                    breakpoint_dict[seg[0]].append(seg[1])
                found = False
                for bkp in breakpoint_dict[seg[0]]:
                    if abs(bkp - seg[2]) < delta:
                        found = True
                        break
                if found is False:
                    breakpoint_dict[seg[0]].append(seg[2])
    for key in breakpoint_dict:
        breakpoint_dict[key].sort()
    splitSimHgtList = []
    for sim_seq in simHgtList:
        tmp_sim_seq = []
        for sim_seg in sim_seq:
            candidate_split_bkp = []
            for bkp in breakpoint_dict[sim_seg[0]]:
                if bkp - sim_seg[1] > delta and sim_seg[2] - bkp > delta:
                    candidate_split_bkp.append(bkp)
            if len(candidate_split_bkp) == 0:
                tmp_sim_seq.append(sim_seg)
            else:
                candidate_split_bkp.insert(0, sim_seg[1])
                candidate_split_bkp.append(sim_seg[2])
                candidate_split_bkp.sort()
                for i in range(0, len(candidate_split_bkp) - 1):
                    tmp_seg = [sim_seg[0], candidate_split_bkp[i], candidate_split_bkp[i + 1], sim_seg[-1]]
                    tmp_sim_seq.append(tmp_seg)
        splitSimHgtList.append(tmp_sim_seq)
    splitTrueHGT = []
    for true_seg in trueHGT:
        candidate_split_bkp = []
        for bkp in breakpoint_dict[true_seg[0]]:
            if bkp - true_seg[1] > delta and true_seg[2] - bkp > delta:
                candidate_split_bkp.append(bkp)
        if len(candidate_split_bkp) == 0:
            splitTrueHGT.append(true_seg)
        else:
            candidate_split_bkp.insert(0, true_seg[1])
            candidate_split_bkp.append(true_seg[2])
            candidate_split_bkp.sort()
            for i in range(0, len(candidate_split_bkp) - 1):
                tmp_seg = [true_seg[0], candidate_split_bkp[i], candidate_split_bkp[i + 1], true_seg[-1]]
                splitTrueHGT.append(tmp_seg)
    return splitTrueHGT, splitSimHgtList


def calc_score(matrix, x, y, seq1, seq2):
    similarity  =   match if seq1[x-1][0] == seq2[y-1][0] and abs(seq1[x-1][1] - seq2[y-1][1]) < delta and abs(seq1[x-1][2] - seq2[y-1][2]) < delta else mismatch
    diag_score  =   matrix[x - 1][y - 1] + similarity
    up_score    =   matrix[x - 1][y] + gap
    left_score  =   matrix[x][y - 1] + gap
    return max(0, diag_score, up_score, left_score)


def create_score_matrix(rows, cols, seq1, seq2):
    score_matrix = [[0 for col in range(cols)] for row in range(rows)]
    max_score = 0
    max_pos   = None
    for i in range(1, rows):
        for j in range(1,cols):
            score = calc_score(score_matrix, i, j, seq1, seq2)
            if score > max_score:
                max_score = score
                max_pos = (i,j)
            score_matrix[i][j] = score
    return score_matrix, max_pos, max_score



def findDetectedHGTSegment(receptor_name, reconstructed_hgt_seq):
    ls = reconstructed_hgt_seq
    hgt_seg_list = []
    i = 0
    while i < len(ls):
        if ls[i][0] != receptor_name:
            hgt_seg_list.append(ls[i][0])
        i += 1
    return hgt_seg_list





def outputReconstructResult(outfile, true_hgt_seg_dict, detected_HGT_seg_dict):
    fo = open(outfile, 'a')
    reconstruct_count = 0
    for receptor_name in detected_HGT_seg_dict:
        reconstruct_count += len(detected_HGT_seg_dict[receptor_name])
    fo.write(str(order)+','+str(hgt_count)+','+str(reconstruct_count)+'\n')
    fo.close()


parser = argparse.ArgumentParser(description="Get HGT reference", add_help=False, usage="%(prog)s [-h] -r genome_dir -s sample_name.txt", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
required = parser.add_argument_group("required arguments")
optional = parser.add_argument_group("optional arguments")
required.add_argument("-r", type=str, help="<str> Metagenomic reference.", metavar="\b")
required.add_argument("-c", type=str, help="<str> coverage file.", metavar="\b")
required.add_argument("-s", type=str, help="<str> sample name.", metavar="\b")
required.add_argument("-a", type=str, help="<str> accurate breakpoints file.", metavar="\b")
required.add_argument("-o", type=str, default='.', help="<str> path to the directory where result should be stored.", metavar="\b")
optional.add_argument("-h", "--help", action="help")
args = vars(parser.parse_args())


ref_genome = pysam.Fastafile(args["r"])
sample_id = args["s"]
coverage_file = args["c"]
break_point_file = args["a"]

rlen = 100
delta = 20
species_num = 640
#simulation_ref_genome = '/disk2/workspace/lichen/li/big_bkp_sim/simulation_ref_genome.fa'
#ref_genome = pysam.Fastafile(simulation_ref_genome)


match = 2
mismatch = -1
gap = 0



sim_cov_evalue_dict = {}
sim_cov_avg_score_dict = {}
sim_num_harbor_dict = {}

file_name = args["o"] + '/' + sample_id + '_true_hgt.txt'
outfile = args["o"] + '/' + sample_id + '_detected_hgt_rate.txt' 
reconstructed_hgt_file = args["o"] + '/' + sample_id+'_reconstructed.txt'
true_hgt_dict, true_hgt_seg_dict = readTruthSeq(file_name)
detected_HGT_seg_dict = {}
cov_evalue_dict = {}
cov_avg_score = {}
num_harbor_dict = {}
all_seq_dict = {}
replication_list = []
seq_dict_list = []
candidate_reference_dict_list = []
ref_coverage_dict = {}
seg_coverage_dict = {}
remained_coverage_dict = {}
ref_name_list = []
fi = open(break_point_file, "r")
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    buf = re.split(r'[:;,\s]\s*', buf)
    from_ref = buf[0]
    to_ref = buf[2]
    if from_ref not in ref_name_list:
        ref_name_list.append(from_ref)
    if to_ref not in ref_name_list:
        ref_name_list.append(to_ref)
genome_sequence_coverage = getGenomeSeqCoverage(coverage_file, ref_name_list)
sim_seq_length_dict = {}
for ref_name in ref_name_list:
    if ref_name not in sim_seq_length_dict:
        tmp_dict = {ref_name : len(ref_genome.fetch(ref_name, 0, ))}
        sim_seq_length_dict.update(tmp_dict)
junction_edge_list = modify_edge_breakpoint(break_point_file, genome_sequence_coverage)
exclude_seq_name_list = []
all_seq_dict = {}
hgt_seq_dict = {}
order = 0
while 1:
    directed_contact_breakpoint_dict = getDirectedContactBreakpoint(junction_edge_list)
    ref_bkp_dict = get_bkp_dict(directed_contact_breakpoint_dict)
    ConnectRefOfEachRefDict = getConnectRefOfEachRef(ref_bkp_dict, directed_contact_breakpoint_dict)
    for harbor_seq_name in ConnectRefOfEachRefDict:
        if harbor_seq_name not in exclude_seq_name_list:
            exclude_seq_name_list.append(harbor_seq_name)
    candidate_reference_dict = getAssembleReference(ConnectRefOfEachRefDict, order)
    if any(candidate_reference_dict) is False:
        break
    for seq_name in candidate_reference_dict:
        sim_seq_length_dict.update({seq_name : getSeqLength(candidate_reference_dict[seq_name])})
    for seq_name in candidate_reference_dict:
        seq = candidate_reference_dict[seq_name]
        ori_seq = getOriginalComponentOfSeq(seq, all_seq_dict, genome_sequence_coverage)
        all_seq_dict.update({seq_name : ori_seq})
    junction_edge_list, bkp_change_map = updateBreakpoints(all_seq_dict, candidate_reference_dict, directed_contact_breakpoint_dict, ConnectRefOfEachRefDict)
    order += 1
for seq_name in all_seq_dict:
    if seq_name not in exclude_seq_name_list:
        hgt_seq_dict.update({seq_name:all_seq_dict[seq_name]})
all_bkp_dict = getAllBkp(hgt_seq_dict)
seq2harbor_dict, harbor2seq_dict = getHarborNameOfSeq(hgt_seq_dict)
harbor_seg_dict = getHarborSeq(harbor2seq_dict, all_bkp_dict)
hgt_seg_dict = getAllHgtSeg(hgt_seq_dict, seq2harbor_dict)
seg2seq_dict, seg_cov_dict, seq_cov_dict = mapSeg2Seq(all_bkp_dict, hgt_seq_dict, hgt_seg_dict, seq2harbor_dict)
out_junction_edge_dict, back_junction_edge_dict, near_junction_edge_dict = getJunctionEdge(harbor2seq_dict, seq2harbor_dict, harbor_seg_dict, hgt_seq_dict)
junction_seg_cov_dict = getCoverageOfSeq(seg2seq_dict, harbor2seq_dict, harbor_seg_dict, seg_cov_dict, seq_cov_dict, hgt_seg_dict, out_junction_edge_dict, back_junction_edge_dict, near_junction_edge_dict)
junction_graph_dict = getJunctionGraph(harbor2seq_dict, harbor_seg_dict, hgt_seq_dict)
all_possible_hgt_dict, all_possible_hgt_freq_dict = getAllPossibleSeq(junction_graph_dict, junction_seg_cov_dict, harbor_seg_dict, hgt_seq_dict)
for harbor_name in all_possible_hgt_dict:
    for i in range(len(all_possible_hgt_dict[harbor_name])):
        fo = open(reconstructed_hgt_file, 'a')
        fo.write(harbor_name+'_'+str(i)+',frequency:'+str(all_possible_hgt_freq_dict[harbor_name][i])+'\n')
        for seg in all_possible_hgt_dict[harbor_name][i]:
    #splitTrueHGT, splitSimHgtList = split_segment(true_hgt_dict[harbor_name], [all_possible_hgt_dict[harbor_name][max_length_index]])
            fo.write(seg[0]+','+str(seg[1])+','+str(seg[2])+','+seg[3]+'\n')
        fo.close()

