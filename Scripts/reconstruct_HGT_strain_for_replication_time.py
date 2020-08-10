import matplotlib
matplotlib.use('Agg')
import argparse
import sys
import fnmatch
import os
import re
import pysam
import copy
import random
import copy
import numpy as np
from scipy import signal
from scipy import ndimage
from scipy import stats
import lmfit
from itertools import product
import matplotlib.pyplot as plt
import pandas as pd
from gurobipy import *
plt.switch_backend('agg')
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams['pdf.fonttype'] = 42


def gc_skew(seq, window = 1000, slide = 10):
    """
    calculate gc skew and cumulative sum of gc skew over sequence windows
     gc skew = ((G - C) / (G + C)) * window size * genome length    
    """
    # convert to G - C
    replacements = {'G':1, 'C':-1, 'A':0, 'T':0, 'N':0}
    gmc = [] # G - C
    for base in seq:
        try:
            gmc.append(replacements[base])
        except:
            gmc.append(0)
    # convert to G + C
    gpc = [abs(i) for i in gmc] # G + C
    # calculate sliding windows for (G - C) and (G + C)
    weights = np.ones(window)/window
    gmc = [[i, c] for i, c in enumerate(signal.fftconvolve(gmc, weights, 'same').tolist())]
    gpc = [[i, c] for i, c in enumerate(signal.fftconvolve(gpc, weights, 'same').tolist())]
    # calculate gc skew and cummulative gc skew sum
    skew = [[], []] # x and y for gc skew
    c_skew = [[], []] # x and y for gc skew cummulative sums
    cs = 0 # cummulative sum
    # select windows to use based on slide
    for i, m in gmc[0::slide]:
        p = gpc[i][1]
        if p == 0:
            gcs = 0
        else:
            gcs = m/p
        cs += gcs
        skew[0].append(i)
        c_skew[0].append(i)
        skew[1].append(gcs)
        c_skew[1].append(cs)
    #ori, ter = find_ori_ter(c_skew, length)
    return c_skew


def plot_coverage(coverage, cov_label, avg_cov, fit, c_skew, ori, ter, m_filter, ptr, title):
    """
    plot coverage profiles with fit, ori, and ter
    """
    # remove some points for plotting (approx. 1,000 datapoints)
    N = int(len(coverage[0])/1000)
    if N != 0:
        X, Y = [coverage[0][0::N], coverage[1][0::N]]
    else:
        X, Y = coverage
    # plot
    fig, ax1 = plt.subplots()
    # plot coverage data
    ax1.plot(X, Y, label = 'coverage', c = '0.60', marker = 'o', ms = 2)
    ax1.axvline(x = ori[0], c = 'r', \
        label = 'consensus Ori: %s' % (ori[0]), linewidth = 2)
    ax1.axvline(x = ter[0], c = 'b', \
        label = 'consensus Ter: %s' % (ter[0]), linewidth = 2)
    # plot median filter and estimated ori and ter
    if m_filter is not None and m_filter is not False:
        ax1.plot(m_filter[0], m_filter[1], \
            label = 'median filter', c = 'k', linewidth = 3)
    # plot fit
    if fit is not False:
        ax1.plot(coverage[0], fit, label = 'least squares fitting', \
            c = 'm', alpha = 0.75, linewidth = 2)
    # plot cumulative gc skew
    if c_skew is not False:
        ax2 = ax1.twinx()
        ax2.set_ylabel('cumulative GC skew')
        ax2.plot(c_skew[0], c_skew[1], label = 'cumulative GC skew', \
            c = 'g', alpha = 0.75, linewidth = 2)
    # title
    plt.suptitle(title, fontsize = 12)
    if ptr != 'n/a':
        ptr = '%.2f' % (ptr)
    plt.title('ptr: %s     avg. cov: %.2f' % (ptr, avg_cov), fontsize = 10)
    # label axes
    ylab = cov_label
    xlab = 'position on genome (bp)'
    ax1.set_ylabel(ylab)
    ax1.set_xlabel(xlab)
    ax1.set_xlim(min(X), max(X))
    # legend
    ax1.legend(loc = 'upper right', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    plt.legend(loc = 'upper left', \
            bbox_to_anchor=(0.5, -0.1), prop = {'size':10})
    # save
    plt.close()
    return fig


def plot_genomes(sample, plot_name, harbor_name, sub_name):
    """
    plot coverage and fitting data for each genome and sample pair
    """
    # PDF for plost
    if '.pdf' not in plot_name:
        plot_name = '%s.pdf' % (plot_name)
    pdf = PdfPages(plot_name)
    fit = sample['fit']
    c_skew = sample['c_skew']
    # plot coverage
    y = sample['cov']
    x = list(range(0, len(y)))
    fig = plot_coverage([x, y], 'coverage', \
            sample['avg_cov'], False, c_skew, \
            sample['ORI'], sample['TER'], None, \
            sample['ptr'], 'genome: %s sample: %s' % \
            (harbor_name, sub_name))
    pdf.savefig(fig, bbox_inches = 'tight')
    # plot filtered/transformed coverage and fitting
    if sample['filtered'] is not False:
        fig = plot_coverage(sample['filtered'], 'coverage (log2, filtered)', \
                sample['avg_cov'], fit, c_skew, \
                sample['ORI'], sample['TER'], sample['m_filter'], \
                sample['ptr'], 'genome: %s hgt: %s' % \
                (harbor_name, sub_name))
        pdf.savefig(fig, bbox_inches = 'tight')
    # save PDF
    pdf.close()


def log_trans(array):
    """
    log transform elements in array
    - leave 0 as 0
    """
    lt = []
    eps = 1e-50
    for i in array:
        if i < eps:
            lt.append(np.log2(eps))
        else:
            lt.append(np.log2(i))
    return lt


def find_y(X, x, y):
    """
    find y value for item in x closest to X
    """
    return y[sorted([[abs(X-p), i] for i, p in enumerate(x)])[0][1]]


def check_peaks(peaks, length):
    """
    select pair from peaks and troughs that are not too close or
    too far apart and have greatest y distance between one another
    """
    # if ori/ter peaks are too close or too far apart, they are probably wrong
    closest, farthest = int(length * float(0.45)), int(length * float(0.55))
    pairs = []
    for pair in list(product(*peaks)):
        pk, tr = pair # peak and trough
        a = (tr[0] - pk[0]) % length
        b = (pk[0] - tr[0]) % length
        if pk[1] < tr[1]:
            continue
        peak_dist = abs(pk[1] - tr[1]) # distance between values
        if (a <= farthest and a >= closest) or (b <=farthest and b >= closest):
            pairs.append([peak_dist, pk, tr])
    if len(pairs) == 0:
        return False, False
    peak_dist, pk, tr = sorted(pairs, reverse = True)[0]
    return pk, tr


def estimate_pars(pars, window = 999):
    """
    estimate parameters for curve fitting
    """
    genome, sample, xy, length = pars
    if xy is False:
        return False
    x, y = xy
    y_med = [x, median_filter(y)]
    # find indexes of peaks and troughs
    pks = signal.find_peaks_cwt(y, np.arange(100,1000,10000))
    trs = signal.find_peaks_cwt([-i for i in y], np.arange(100,1000,10000))
    # find positions on genome for peaks and troughs
    pks = [[y_med[0][i], y_med[1][i]] for i in pks] 
    trs = [[y_med[0][i], y_med[1][i]] for i in trs]
    # find best pk/tr pair based on greatest distance in coverage 
    # and position on genome
    ori, ter = check_peaks([pks, trs], length)
    x1, x2 = ori[0], ter[0]
    y1, y2 = ori[1], ter[1]
    if genome is not None:
        return genome, sample, (x1, x2, y1, y2, y_med)
    else:
        return x1, x2, y1, y2, y_med


def fit_coverage(pars, window = 999, est_pars = False):
    """
    use coverage data to estimate parameters for
    coverage function
    """
    x, y, length = pars
    # estimate variables using median filter
    if est_pars is True:
        x1_est, x2_est, y1_est, y2_est, y_med = \
            estimate_pars((None, None, (x, y), length), window = window)
    else:
        x1_est = int(length * float(0.25))
        x2_est = int(length * float(0.75))
        y1_est, y2_est = min(y), max(y)
        y_med = [x, median_filter(y)]
    # how close can origin and terminus of 
    # replication be to one another?
    closest, farthest = \
            int(length * float(0.45)), int(length * float(0.55))
    # Parameters
    Pars = lmfit.Parameters()
    ## x1 and x2
    Pars.add('length', value = length, vary = False)
    Pars.add('x2', value = x2_est, min = 0, max = length)
    Pars.add('xdist', value = abs(x1_est - x2_est), min = closest, max = farthest)
    Pars.add('x1', expr = '(x2 + xdist) % length')
    ## y1 and y2
    Pars.add('y1', value = y1_est)
    Pars.add('y2', value = y2_est)
    # fit data to model
    mi = lmfit.minimize(coverage_function, Pars, args = (x,), \
            kws = {'data':y, 'printPs':False}, method = 'leastsq')
    # get fitted values
    fit = coverage_function(mi.params, x)
    return (x1_est, x2_est, y_med, mi.params, fit, mi.redchi)


def coverage_function(pars, X, data = None, printPs = False): 
    """
    piecewise linear function representing 
    coverage profile across genome
    """
    results = []
    x1 = pars['x1'].value
    x2 = pars['x2'].value
    y1, y2 = pars['y1'].value, pars['y2'].value
    if y1 > y2: # y1 ~ ori
        a = float(y2 - y1) / float(x2 - x1)
    else:
        a = float(y1 - y2) / float(x1 - x2)
    if printPs is True:
        print('x1: %s x2: %s y1: %s y2: %s a:%s' \
                % ('{:,}'.format(int(x1)), '{:,}'.format(int(x2)), y1, y2, a))
    for x in X:
        if x <= x1:
            results.append(-1*a*x + y1 + a*x1)
        elif x1 < x < x2:
            results.append(a*x + y1 - a*x1)
        elif x >= x2:
            results.append(-1*a*x + y2 + a*x2)
    if data is None:
        return np.asarray(results)
    return np.asarray([y - data[i] for i, y in enumerate(results)]) # model - data


def ori_from_cov(sample, length, error_threshold = 20000):
    """
    find x, y for ORI and TER based on
    result of fitting
    """
    x = sample[0]
    y = sample[1]
    # find origin and terminus from coverage
    x1_est, x2_est, y_med, pars, fit, chi = \
            fit_coverage([x, y, length], est_pars = False)
    x1, x2 = int(pars['x1'].value), int(pars['x2'].value)
    x1_err, x2_err = pars['x1'].stderr, pars['x2'].stderr
    y1 = find_y(x1, y_med[0], y_med[1])
    y2 = find_y(x2, y_med[0], y_med[1])
    if y1 > y2:
        ORI = (x1, y1, x1_err)
        TER = (x2, y2, x2_err)
    else:
        ORI = (x2, y2, x2_err)
        TER = (x1, y1, x1_err)
    # calculate PTR from coverage at Ori and Ter
    # exclude if x1 or x2 error > error_threshold
    if x1_err > error_threshold or x2_err > error_threshold:
        ptr = False
    elif TER[1] == 0:
        ptr = False
    else:
        ptr = (2**ORI[1])/(2**TER[1])
    m_filter = y_med
    return ptr, ORI, TER, m_filter, fit


def detectReplicatingRef(candidate_reference_dict, seq_dict, sample_id):
    HGT_ReplicatingRefIndexList = []
    non_HGT_ReplicatingRefIndexList = []
    for ref_name in candidate_reference_dict:
        tmp_coverage_list = []
        for position in range(0, sim_seq_length_dict[ref_name]):
            if position in genome_sequence_coverage[ref_name]:
                tmp_coverage_list.append(genome_sequence_coverage[ref_name][position])
        tmp_coverage_array0 = np.asarray(tmp_coverage_list)
        filtered_coverage0 = coverage_windows(tmp_coverage_array0)
        if filtered_coverage0 is False:
            continue
        avg_cov0 = np.average(filtered_coverage0[1])
        filtered_coverage0[1] = log_trans(filtered_coverage0[1])
        seq = ref_genome.fetch(ref_name, 0, sim_seq_length_dict[ref_name])
        ptr0, ORI0, TER0, m_filter0, fit0 = ori_from_cov(filtered_coverage0, sim_seq_length_dict[ref_name], error_threshold = 20000)
        c_skew0 = gc_skew(seq, window = 1000, slide = 10)
        flag = False
        if ptr0 is False:
            continue
        for i in range(0, len(candidate_reference_dict[ref_name])):
            candidate_reference = candidate_reference_dict[ref_name][i]
            candidate_reference_length = 0
            for segment in candidate_reference:
                candidate_reference_length += abs(segment[1] - segment[2])
            tmp_coverage_list = []
            seq = seq_dict[ref_name][i]
            for segment in candidate_reference:
                seg_name = segment[0]
                seg_left_pos = min(segment[1], segment[2])
                seg_right_pos = max(segment[1], segment[2])
                sign = segment[3]
                if sign is '+':
                    for position in range(seg_left_pos, seg_right_pos):
                        if position in genome_sequence_coverage[seg_name]:
                            tmp_coverage_list.append(genome_sequence_coverage[seg_name][position])
                else:
                    for position in range(seg_right_pos, seg_left_pos, -1):
                        if position in genome_sequence_coverage[seg_name]:
                            tmp_coverage_list.append(genome_sequence_coverage[seg_name][position])
            tmp_coverage_array = np.asarray(tmp_coverage_list)
            filtered_coverage = coverage_windows(tmp_coverage_array)
            if filtered_coverage is False:
                continue
            avg_cov = np.average(filtered_coverage[1])
            filtered_coverage[1] = log_trans(filtered_coverage[1])
            ptr, ORI, TER, m_filter, fit = ori_from_cov(filtered_coverage, candidate_reference_length, error_threshold = 20000)
            c_skew = gc_skew(seq, window = 1000, slide = 10)
            if ptr is False:
                continue
            elif max(ptr0, ptr) / min(ptr0, ptr) > 2:
                continue
            sample = {'c_skew' : c_skew, 'fit' : fit, 'cov' : tmp_coverage_array, 'ptr' : ptr, 'avg_cov' : avg_cov, 'ORI' : ORI, 'TER' : TER, 'filtered' : filtered_coverage, 'm_filter' : m_filter}
            plot_name = sample_id + '_' + ref_name + '_hgt_' + str(i)
            plot_genomes(sample, plot_name, ref_name, str(i))
            HGT_ReplicatingRefIndexList.append([ref_name, i, ptr])
            flag = True
        if flag is False:
            sample = {'c_skew' : c_skew0, 'fit' : fit0, 'cov' : tmp_coverage_array0, 'ptr' : ptr0, 'avg_cov' : avg_cov0, 'ORI' : ORI0, 'TER' : TER0, 'filtered' : filtered_coverage0, 'm_filter' : m_filter0}
            plot_name = sample_id + '_' + ref_name
            plot_genomes(sample, plot_name, ref_name, 'non_hgt')
            non_HGT_ReplicatingRefIndexList.append([ref_name, ptr0])
    return HGT_ReplicatingRefIndexList, non_HGT_ReplicatingRefIndexList



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
        if len(buf) > 7:
            from_ref = buf[0]
            from_pos = int(buf[1])
            to_ref = buf[4]
            to_pos = int(buf[5])
        else:
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
        from_side = buf[-3]
        if len(buf) > 7:
            to_ref = buf[4]
            to_pos = int(buf[5])
        else:
            to_ref = buf[2]
            to_pos = int(buf[3])
        to_pos = change_map_dict[to_ref][to_pos]
        to_side = buf[-2]
        if buf[-1] == 'False':
            reverse = False
        else:
            reverse = True
        if from_ref not in genome_sequence_coverage or to_ref not in genome_sequence_coverage:
            continue
        edge = [from_ref, from_pos, from_side, to_ref, to_pos, to_side, reverse]
        junction_edge_list.append(edge)
    return junction_edge_list


def getDirectedContactBreakpoint(junction_edge_list, sim_seq_length_dict):
    directed_contact_breakpoint_dict = {}
    for edge in junction_edge_list:
        ref1 = edge[0]
        ref2 = edge[3]
        bkp1 = edge[1]
        bkp2 = edge[4]
        bkp1_side = edge[2]
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
                        if breakpoint[0] == to_pos and breakpoint[3] == reverse and breakpoint[1] == from_side and breakpoint[2] == to_side:
                            found = True
                            break
                    if found is False:
                        directed_contact_breakpoint_dict[from_ref][from_pos][to_ref].append([to_pos, from_side, to_side, reverse])
        if from_ref not in sim_seq_length_dict:
            tmp_dict = {from_ref : len(ref_genome.fetch(from_ref, 0, ))}
            sim_seq_length_dict.update(tmp_dict)
        if to_ref not in sim_seq_length_dict:
            tmp_dict = {to_ref : len(ref_genome.fetch(to_ref, 0, ))}
            sim_seq_length_dict.update(tmp_dict)
    return directed_contact_breakpoint_dict, sim_seq_length_dict



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
        min_bkp = min(candidate_bkp_list)
        if ref_name in directed_contact_breakpoint_dict:
            if min_bkp in directed_contact_breakpoint_dict[ref_name]:
                if abs(0 - min(ref_bkp_dict[ref_name])) > 150:
                    ref_bkp_dict[ref_name].append(0)
        ref_bkp_dict[ref_name].sort()
    return ref_bkp_dict



def prepareGeneLibrary(directed_contact_breakpoint_dict, sim_seq_length_dict, ref_bkp_dict):
    ref_segment_dict = {}
    for ref_name in ref_bkp_dict:
        if len(ref_bkp_dict[ref_name]) <= 1:
            continue
        ls = []
        count = 1
        end_pos = sim_seq_length_dict[ref_name]
        for i in range(0, len(ref_bkp_dict[ref_name]) - 1):
            left_pos = ref_bkp_dict[ref_name][i]
            right_pos = ref_bkp_dict[ref_name][i + 1]
            found = False
            if ref_name in directed_contact_breakpoint_dict:
                if left_pos in directed_contact_breakpoint_dict[ref_name]:
                    for to_ref in directed_contact_breakpoint_dict[ref_name][left_pos]:
                        if to_ref not in ref_bkp_dict:
                            continue
                        for bkp_info in directed_contact_breakpoint_dict[ref_name][left_pos][to_ref]:
                            if bkp_info[1] == 'left':
                                found = True
                                break
                        if found is True:
                            break
            if found is False:
                for from_ref in directed_contact_breakpoint_dict:
                    if from_ref not in ref_bkp_dict:
                        continue
                    for from_pos in directed_contact_breakpoint_dict[from_ref]:
                        if ref_name in directed_contact_breakpoint_dict[from_ref][from_pos]:
                            for bkp_info in directed_contact_breakpoint_dict[from_ref][from_pos][ref_name]:
                                if bkp_info[0] == left_pos and bkp_info[2] == 'left':
                                    found = True
                                    break
                        if found is True:
                            break
                    if found is True:
                        break
            if found is False:
                if ref_name in directed_contact_breakpoint_dict:
                    if right_pos in directed_contact_breakpoint_dict[ref_name]:
                        for to_ref in directed_contact_breakpoint_dict[ref_name][right_pos]:
                            if to_ref not in ref_bkp_dict:
                                continue
                            for bkp_info in directed_contact_breakpoint_dict[ref_name][right_pos][to_ref]: 
                                if bkp_info[1] == 'right':
                                    found = True
                                    break
                            if found is True:
                                break
                if found is False:
                    for from_ref in directed_contact_breakpoint_dict:
                        if from_ref not in ref_bkp_dict:
                            continue
                        for from_pos in directed_contact_breakpoint_dict[from_ref]:
                            if ref_name in directed_contact_breakpoint_dict[from_ref][from_pos]:
                                for bkp_info in directed_contact_breakpoint_dict[from_ref][from_pos][ref_name]:
                                    if bkp_info[0] == right_pos and bkp_info[2] == 'right':
                                        found = True
                                        break
                            if found is True:
                                break
                        if found is True:
                            break
            if found is True:
                tmp = [count, left_pos, right_pos, False]
                ls.append(tmp)
                count += 1
        left_pos = ref_bkp_dict[ref_name][-1]
        right_pos = end_pos
        if right_pos - left_pos > delta:
            found = False
            if ref_name in directed_contact_breakpoint_dict:
                if left_pos in directed_contact_breakpoint_dict[ref_name]:
                    for to_ref in directed_contact_breakpoint_dict[ref_name][left_pos]:
                        if to_ref not in ref_bkp_dict:
                            continue
                        for bkp_info in directed_contact_breakpoint_dict[ref_name][left_pos][to_ref]:
                            if bkp_info[1] == 'left':
                                found = True
                                break
                        if found is True:
                            break
            if found is False:
                for from_ref in directed_contact_breakpoint_dict:
                    if from_ref not in ref_bkp_dict:
                        continue
                    for from_pos in directed_contact_breakpoint_dict[from_ref]:
                        if ref_name in directed_contact_breakpoint_dict[from_ref][from_pos]:
                            for bkp_info in directed_contact_breakpoint_dict[from_ref][from_pos][ref_name]:
                                if bkp_info[0] == left_pos and bkp_info[2] == 'left':
                                    found = True
                                    break
                        if found is True:
                            break
                    if found is True:
                        break
            if found is True:
                tmp = [count, left_pos, right_pos, True]
                ls.append(tmp)
        if len(ls) > 0:
            if ls[0][1] != 0:
                tmp_seg = [1, 0, ls[0][1], False]
                for j in range(0, len(ls)):
                    ls[j][0] += 1
                ls.insert(0, tmp_seg)
            tmp_dict = {ref_name : ls}
            ref_segment_dict.update(tmp_dict)
    return ref_segment_dict



def getDiffRef(currentVisitedRefList):
    visited_ref = []
    for junction_edge in currentVisitedRefList:
        if junction_edge[0] not in visited_ref:
            visited_ref.append(junction_edge[0])
        if junction_edge[1] not in visited_ref:
            visited_ref.append(junction_edge[1])
    return len(visited_ref)


def BFS_modify(from_ref, directed_contact_breakpoint_dict, ref_num_threshold, currentVisitedRefList, target_ref, connectedRefList, ref_segment_dict):
    num_diff_ref = getDiffRef(currentVisitedRefList)
    if len(currentVisitedRefList) > 1:
        return connectedRefList
    if num_diff_ref < ref_num_threshold:
        for from_pos in directed_contact_breakpoint_dict[from_ref]:
            for to_ref in directed_contact_breakpoint_dict[from_ref][from_pos]:
                if len(currentVisitedRefList) == 1:
                    if currentVisitedRefList[-1][0] != to_ref:
                        continue
                if to_ref not in ref_segment_dict:
                    continue
                for segment in directed_contact_breakpoint_dict[from_ref][from_pos][to_ref]:
                    from_side = segment[-3]
                    to_side = segment[-2]
                    reverse = segment[-1]
                    to_pos = segment[0]
                    flag = True
                    if len(currentVisitedRefList) == 1:
                        previous_to_side = currentVisitedRefList[-1][-2]
                        previous_reverse = currentVisitedRefList[-1][-1]
                        if previous_to_side == 'right':
                            if from_side != 'left':
                                continue
                        else:
                            if from_side != 'right':
                                continue
                        if currentVisitedRefList[-1][0] == to_ref:
                            if reverse != previous_reverse:
                                continue
                        if to_ref == currentVisitedRefList[-1][0]:
                            if reverse is False:
                                if to_pos < currentVisitedRefList[-1][-5]:
                                    flag = False
                                    break
                            insert_length = abs(currentVisitedRefList[-1][-4] - from_pos)
                            bkp_length = abs(currentVisitedRefList[-1][-5] - to_pos)
                            if bkp_length > 300:
                                flag = False
                                break
                            from_coverage = getRefCoverage(from_ref)
                            from_seg_coverage = getSegCoverage(from_ref, min(currentVisitedRefList[-1][-4], from_pos), max(currentVisitedRefList[-1][-4], from_pos))
                            to_coverage = getRefCoverage(to_ref)
                            if insert_length < 300:
                                if from_seg_coverage >= (from_coverage + to_coverage) * 0.7:
                                    if reverse is False:
                                        if currentVisitedRefList[-1][-4] > from_pos:
                                            flag = False
                                            break
                                    harbor_ref_name = to_ref
                                    harbor_ref_length = sim_seq_length_dict[harbor_ref_name]
                                    if harbor_ref_length + insert_length > 20000000:
                                        flag = False
                                        break
                                else:
                                    flag = False
                                    if min(from_coverage, to_coverage) / max(from_coverage, to_coverage) > 0.8:
                                        if (min(insert_length, bkp_length) > 0 and min(insert_length, bkp_length) / max(insert_length, bkp_length) > 0.7) or insert_length <= 150 or bkp_length <= 150:
                                            flag = True
                                            if reverse is False:
                                                if currentVisitedRefList[-1][-4] < from_pos:
                                                    flag = False
                                                    break
                                            short_ref_length = min(sim_seq_length_dict[from_ref], sim_seq_length_dict[to_ref])
                                            long_ref_length = max(sim_seq_length_dict[from_ref], sim_seq_length_dict[to_ref])
                                            if short_ref_length + long_ref_length > 20000000:
                                                flag = False
                                                break
                                            if short_ref_length / long_ref_length > 0.5:
                                                flag = False
                                                break
                            else:
                                if from_seg_coverage < (from_coverage + to_coverage) * 0.7:
                                    flag = False
                                    break
                                if reverse is False:
                                    if currentVisitedRefList[-1][-4] > from_pos:
                                        flag = False
                                        break
                                harbor_ref_name = to_ref
                                harbor_ref_length = sim_seq_length_dict[harbor_ref_name]
                                if harbor_ref_length + insert_length > 20000000:
                                    flag = False
                                    break
                    if flag is False:
                        continue
                    if to_ref == target_ref:
                        currentVisitedRefList.append([from_ref, to_ref, from_pos, to_pos, from_side, to_side, reverse])#from_ref, to_ref, from_pos, tp_pos
                        tmp_ls = copy.deepcopy(currentVisitedRefList)
                        connectedRefList.append(tmp_ls)
                        currentVisitedRefList.pop()
                        return connectedRefList
                    else:
                        currentVisitedRefList.append([from_ref, to_ref, from_pos, to_pos, from_side, to_side, reverse])
                        connectedRefList = BFS_modify(to_ref, directed_contact_breakpoint_dict, ref_num_threshold, currentVisitedRefList, target_ref, connectedRefList, ref_segment_dict)
                        currentVisitedRefList.pop()
    return connectedRefList



def getConnectRefOfEachRef(ref_segment_dict, directed_contact_breakpoint_dict):
    ConnectRefOfEachRefDict = {}
    ref_num_threshold = 3
    count = 0
    for from_ref in ref_segment_dict:
        if ref_segment_dict[from_ref][0][1] != 0:
            continue
        count += 1
    count = 0
    for from_ref in ref_segment_dict:
        if ref_segment_dict[from_ref][0][1] != 0:
            continue
        if from_ref in directed_contact_breakpoint_dict:
            count += 1
            connectedRefList = []
            currentVisitedRefList = []
            connectedRefList = BFS_modify(from_ref, directed_contact_breakpoint_dict, ref_num_threshold, currentVisitedRefList, from_ref, connectedRefList, ref_segment_dict)
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



def getTargetSeg(ref_name, out_back_info, candidate_reference, start_pos, seq, hgt_segment_list):
    interact_ref_name = out_back_info['interact_ref_name']
    current_out_pos = out_back_info['out_pos']
    current_back_pos = out_back_info['back_pos']
    target_pos = out_back_info['target_pos']
    source_pos = out_back_info['source_pos']
    #if abs(source_pos - target_pos) > 10000:
    if out_back_info['out_side'] == 'right':#from ref is '+'
        if out_back_info['out_reverse'] == False:
            if source_pos > target_pos:
                from_coverage_1 = getSegCoverage(ref_name, max(0, min(current_out_pos, current_back_pos) - abs(source_pos - target_pos)), min(current_out_pos, current_back_pos))
                from_coverage_2 = getSegCoverage(ref_name, min(current_out_pos, current_back_pos), min(sim_seq_length_dict[ref_name], min(current_out_pos, current_back_pos) + abs(source_pos - target_pos)))
                from_coverage = (from_coverage_1 + from_coverage_2) / 2
                to_coverage, to_left_pos, to_right_pos = consumeAndGetSegCoverage(interact_ref_name, target_pos, source_pos, from_coverage)
                near_to_coverage_1 = getSegCoverage(interact_ref_name, max(0, target_pos - abs(source_pos - target_pos)), target_pos)
                near_to_coverage_2 = getSegCoverage(interact_ref_name, source_pos, min(source_pos + abs(source_pos - target_pos), sim_seq_length_dict[interact_ref_name]))
                tmp = (near_to_coverage_1 + near_to_coverage_2) / 2
                to_coverage -= tmp
            else:
                from_coverage = getRefCoverage(ref_name)
                to_coverage = consumeAndGetRefCoverage(interact_ref_name, from_coverage)
        else:
            if source_pos < target_pos:
                from_coverage_1 = getSegCoverage(ref_name, max(0, min(current_out_pos, current_back_pos) - abs(source_pos - target_pos)), min(current_out_pos, current_back_pos))
                from_coverage_2 = getSegCoverage(ref_name, min(current_out_pos, current_back_pos), min(sim_seq_length_dict[ref_name], min(current_out_pos, current_back_pos) + abs(source_pos - target_pos)))
                from_coverage = (from_coverage_1 + from_coverage_2) / 2
                to_coverage, to_left_pos, to_right_pos = consumeAndGetSegCoverage(interact_ref_name, source_pos, target_pos, from_coverage)
                near_to_coverage_1 = getSegCoverage(interact_ref_name, max(0, source_pos - abs(source_pos - target_pos)), source_pos)
                near_to_coverage_2 = getSegCoverage(interact_ref_name, target_pos, min(target_pos + abs(source_pos - target_pos), sim_seq_length_dict[interact_ref_name]))
                tmp = (near_to_coverage_1 + near_to_coverage_2) / 2
                to_coverage -= tmp
            else:
                from_coverage = getRefCoverage(ref_name)
                to_coverage = consumeAndGetRefCoverage(interact_ref_name, from_coverage)
    else:
        if out_back_info['out_reverse'] == True:
            if source_pos > target_pos:
                from_coverage_1 = getSegCoverage(ref_name, max(0, min(current_out_pos, current_back_pos) - abs(source_pos - target_pos)), min(current_out_pos, current_back_pos))
                from_coverage_2 = getSegCoverage(ref_name, min(current_out_pos, current_back_pos), min(sim_seq_length_dict[ref_name], min(current_out_pos, current_back_pos) + abs(source_pos - target_pos)))
                from_coverage = (from_coverage_1 + from_coverage_2) / 2
                to_coverage, to_left_pos, to_right_pos = consumeAndGetSegCoverage(interact_ref_name, target_pos, source_pos, from_coverage)
                near_to_coverage_1 = getSegCoverage(interact_ref_name, max(0, target_pos - abs(source_pos - target_pos)), target_pos)
                near_to_coverage_2 = getSegCoverage(interact_ref_name, source_pos, min(source_pos + abs(source_pos - target_pos), sim_seq_length_dict[interact_ref_name]))
                tmp = (near_to_coverage_1 + near_to_coverage_2) / 2
                to_coverage -= tmp
            else:
                from_coverage = getRefCoverage(ref_name)
                to_coverage = consumeAndGetRefCoverage(interact_ref_name, from_coverage)
    if to_coverage <= 0 or to_coverage / from_coverage < 0.5:
        if min(current_out_pos, current_back_pos) > start_pos:
            seq += ref_genome.fetch(ref_name, start_pos, min(current_out_pos, current_back_pos))
            candidate_reference.append([ref_name, start_pos, min(current_out_pos, current_back_pos), '+'])
            current_max_pos = min(current_out_pos, current_back_pos)
        else:
            current_max_pos = min(current_out_pos, current_back_pos)
        if current_out_pos == current_back_pos:
            current_max_pos += 1
        if out_back_info['out_side'] == 'right':#from ref is '+'
            if out_back_info['out_reverse'] == False:
                if source_pos > target_pos:
                    remained_seg_coverage_dict[interact_ref_name][to_left_pos][to_right_pos] += from_coverage
                else:
                    remained_coverage_dict[interact_ref_name] += from_coverage
            else:
                if source_pos < target_pos:
                    remained_seg_coverage_dict[interact_ref_name][to_left_pos][to_right_pos] += from_coverage
                else:
                    remained_coverage_dict[interact_ref_name] += from_coverage
        else:
            if out_back_info['out_reverse'] == True:
                if source_pos > target_pos:
                    remained_seg_coverage_dict[interact_ref_name][to_left_pos][to_right_pos] += from_coverage
                else:
                    remained_coverage_dict[interact_ref_name] += from_coverage
        return candidate_reference, seq, current_max_pos
    if start_pos < min(current_out_pos, current_back_pos):
        seq += ref_genome.fetch(ref_name, start_pos, min(current_out_pos, current_back_pos))
        candidate_reference.append([ref_name, start_pos, min(current_out_pos, current_back_pos), '+'])
    current_max_pos = max(current_out_pos, current_back_pos)
    if out_back_info['out_side'] == 'right':#from ref is '+'
        if out_back_info['out_reverse'] == False:
            if source_pos > target_pos:
                if [interact_ref_name, target_pos, source_pos, '+'] not in hgt_segment_list:
                    candidate_reference.append([interact_ref_name, target_pos, source_pos, '+'])
                    seq += ref_genome.fetch(interact_ref_name, target_pos, source_pos)
                    hgt_segment_list.append([interact_ref_name, target_pos, source_pos, '+'])
            else:
                end_pos = sim_seq_length_dict[interact_ref_name]
                candidate_reference.append([interact_ref_name, target_pos, end_pos, '+'])
                seq += ref_genome.fetch(interact_ref_name, target_pos, end_pos)
                candidate_reference.append([interact_ref_name, 0, source_pos, '+'])
                seq += ref_genome.fetch(interact_ref_name, 0, source_pos)
        else:
            if source_pos < target_pos:
                if [interact_ref_name, source_pos, target_pos, '-'] not in hgt_segment_list:
                    candidate_reference.append([interact_ref_name, source_pos, target_pos, '-'])
                    s = ref_genome.fetch(interact_ref_name, source_pos, target_pos)
                    s = get_mate_seq(s)
                    seq += s
                    hgt_segment_list.append([interact_ref_name, source_pos, target_pos, '-'])
            else:
                end_pos = sim_seq_length_dict[interact_ref_name]
                candidate_reference.append([interact_ref_name, 0, target_pos, '-'])
                s = ref_genome.fetch(interact_ref_name, 0, target_pos)
                s = get_mate_seq(s)
                seq += s
                candidate_reference.append([interact_ref_name, source_pos, end_pos, '-'])
                s = ref_genome.fetch(interact_ref_name, source_pos, end_pos)
                s = get_mate_seq(s)
                seq += s
    else:#from_ref is '-'
        if out_back_info['out_reverse'] == True:
            if source_pos > target_pos:
                if [interact_ref_name, target_pos, source_pos, '-'] not in hgt_segment_list:
                    candidate_reference.append([interact_ref_name, target_pos, source_pos, '-'])
                    s = ref_genome.fetch(interact_ref_name, target_pos, source_pos)
                    s = get_mate_seq(s)
                    seq += s
                    hgt_segment_list.append([interact_ref_name, target_pos, source_pos, '-'])
            else:
                end_pos = sim_seq_length_dict[interact_ref_name]
                s = ref_genome.fetch(interact_ref_name, target_pos, end_pos)
                s1 = get_mate_seq(s)
                s = ref_genome.fetch(interact_ref_name, 0, source_pos)
                s2 = get_mate_seq(s)
                seq += s2
                seq += s1
                candidate_reference.append([interact_ref_name, 0, source_pos, '-'])
                candidate_reference.append([interact_ref_name, target_pos, end_pos, '-'])
    if current_out_pos == current_back_pos:
        current_max_pos += 1
    return candidate_reference, seq, current_max_pos


def getRefNum(candidate_reference):
    ref_name_list = []
    for segment in candidate_reference:
        if segment[0] not in ref_name_list:
            ref_name_list.append(segment[0])
    return len(ref_name_list)


def getAssembleReference(ConnectRefOfEachRefDict):
    candidate_reference_dict = {}
    seq_dict = {}
    for ref_name in ConnectRefOfEachRefDict:
        edge_path_list = ConnectRefOfEachRefDict[ref_name]
        out_back_info_list = []
        for edge_path in edge_path_list:
            if edge_path[0][0] == ref_name and edge_path[1][1] == ref_name:
                out_pos = edge_path[0][2]
                out_side = edge_path[0][-3]
                out_reverse = edge_path[0][-1]
                interact_ref_name = edge_path[0][1]
                target_pos = edge_path[0][3]
                back_pos = edge_path[1][3]
                back_side = edge_path[1][-2]
                back_reverse = edge_path[1][-1]
                source_pos = edge_path[1][2]
                if source_pos > sim_seq_length_dict[interact_ref_name]:
                    source_pos = sim_seq_length_dict[interact_ref_name]
                if target_pos > sim_seq_length_dict[interact_ref_name]:
                    target_pos = sim_seq_length_dict[interact_ref_name]
                if out_pos > sim_seq_length_dict[ref_name]:
                    out_pos = sim_seq_length_dict[ref_name]
                if back_pos > sim_seq_length_dict[ref_name]:
                    back_pos = sim_seq_length_dict[ref_name]
                out_back_info = {'out_pos' : edge_path[0][2], 'out_side' : edge_path[0][-3], 'out_reverse' : edge_path[0][-1], 'interact_ref_name' : edge_path[0][1], 'target_pos' : edge_path[0][3], 'back_pos' : edge_path[1][3], 'back_side' : edge_path[1][-2], 'back_reverse' : edge_path[1][-1], 'source_pos' : edge_path[1][2], 'visited' : False}
                #out_back_info = [out_pos, out_side, out_reverse, target_pos, back_pos, back_side, back_reverse, source_pos, interact_ref_name, False]
                if out_back_info not in out_back_info_list:
                    out_back_info_list.append(out_back_info)
        out_back_info_list.sort(key=lambda ele:ele['out_pos'])
        #out_back_info_list.sort()
        candidate_reference_list = []
        hgt_segment_list = []
        seq_list = []
        for i in range(0, len(out_back_info_list)):
            out_back_info = out_back_info_list[i]
            if out_back_info['visited'] == False:
                candidate_reference = []
                seq = ''
                if out_back_info['out_side'] == 'left' and out_back_info['out_pos'] < out_back_info['back_pos']:
                    continue
                if out_back_info['out_side'] == 'right' and out_back_info['out_pos'] > out_back_info['back_pos']:
                    continue
                candidate_reference, seq, current_max_pos = getTargetSeg(ref_name, out_back_info, candidate_reference, 0, seq, hgt_segment_list)
                out_back_info['visited'] = True
                if i != len(out_back_info_list) - 1:
                    for j in range(i + 1, len(out_back_info_list)):
                        if out_back_info_list[j]['out_pos'] < current_max_pos + 150 or out_back_info_list[j]['back_pos'] < current_max_pos + 150:
                            continue
                        out_back_info = out_back_info_list[j]
                        if out_back_info['visited'] == True:
                            continue
                        if out_back_info['out_side'] == 'left' and out_back_info['out_pos'] < out_back_info['back_pos']:
                            continue
                        if out_back_info['out_side'] == 'right' and out_back_info['out_pos'] > out_back_info['back_pos']:
                            continue
                        candidate_reference, seq, current_max_pos = getTargetSeg(ref_name, out_back_info, candidate_reference, current_max_pos, seq, hgt_segment_list)
                        out_back_info['visited'] = True
                        if len(seq) > 20000000:
                            break
                if len(candidate_reference) != 0:
                    if current_max_pos < sim_seq_length_dict[ref_name]:
                        candidate_reference.append([ref_name, current_max_pos, sim_seq_length_dict[ref_name], '+'])
                        seq += ref_genome.fetch(ref_name, current_max_pos, sim_seq_length_dict[ref_name])
                        if candidate_reference not in candidate_reference_list and len(seq) < 20000000:
                            RefNum = getRefNum(candidate_reference)
                            if RefNum > 1:
                                candidate_reference_list.append(candidate_reference)
                                seq_list.append(seq)
        if len(candidate_reference_list) > 0:
            tmp_dict = {ref_name : candidate_reference_list}
            candidate_reference_dict.update(tmp_dict)
            tmp_dict = {ref_name : seq_list}
            seq_dict.update(tmp_dict)
    return candidate_reference_dict, seq_dict


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
    remained_coverage_dict.update({ref_name : ref_coverage - CostRefCoverage})
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


def consumeAndGetSegCoverage(ref_name, left_pos, right_pos, CostSegCoverage):
    if ref_name not in remained_seg_coverage_dict:
        seg_coverage = getSegCoverage(ref_name, left_pos, right_pos)
        remained_seg_coverage_dict.update({ref_name : {left_pos : {right_pos : seg_coverage - CostSegCoverage}}})
        return seg_coverage, left_pos, right_pos
    else:
        for candidate_left_pos in remained_seg_coverage_dict[ref_name]:
            if abs(left_pos - candidate_left_pos) < 150 and candidate_left_pos < right_pos:
                for candidate_right_pos in remained_seg_coverage_dict[ref_name][candidate_left_pos]:
                    if abs(right_pos - candidate_right_pos) < 150:
                        seg_coverage = remained_seg_coverage_dict[ref_name][candidate_left_pos][candidate_right_pos]
                        remained_seg_coverage_dict[ref_name][candidate_left_pos].update({candidate_right_pos : seg_coverage - CostSegCoverage})
                        return seg_coverage, candidate_left_pos, candidate_right_pos
                seg_coverage = getSegCoverage(ref_name, candidate_left_pos, right_pos)
                remained_seg_coverage_dict[ref_name][candidate_left_pos].update({right_pos : seg_coverage - CostSegCoverage})
                return seg_coverage, candidate_left_pos, right_pos
        seg_coverage = getSegCoverage(ref_name, left_pos, right_pos)
        remained_seg_coverage_dict[ref_name].update({left_pos : {right_pos : seg_coverage - CostSegCoverage}})
        return seg_coverage, left_pos, right_pos
        
        




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
    #sys.exit(main())


delta = 20


replication_list = []
seq_dict_list = []
candidate_reference_dict_list = []

ref_coverage_dict = {}
seg_coverage_dict = {}
remained_coverage_dict = {}
remained_seg_coverage_dict = {}
coverage_file = args["c"]
break_point_file = args["a"]
ref_name_list = []
fi = open(break_point_file, "r")
while 1:
    buf = fi.readline()
    buf = buf.strip('\n')
    if not buf:
        break
    buf = re.split(r'[:;,\s]\s*', buf)
    from_ref = buf[0]
    to_ref = buf[4]
    if from_ref not in ref_name_list:
        ref_name_list.append(from_ref)
    if to_ref not in ref_name_list:
        ref_name_list.append(to_ref)
genome_sequence_coverage = getGenomeSeqCoverage(coverage_file, ref_name_list)
junction_edge_list = modify_edge_breakpoint(break_point_file, genome_sequence_coverage)
sim_seq_length_dict = {}
directed_contact_breakpoint_dict, sim_seq_length_dict = getDirectedContactBreakpoint(junction_edge_list, sim_seq_length_dict)
ref_bkp_dict = get_bkp_dict(directed_contact_breakpoint_dict)
ref_segment_dict = prepareGeneLibrary(directed_contact_breakpoint_dict, sim_seq_length_dict, ref_bkp_dict)
ConnectRefOfEachRefDict = getConnectRefOfEachRef(ref_segment_dict, directed_contact_breakpoint_dict)
#ConnectRefOfEachRefDict_list.append(ConnectRefOfEachRefDict)
candidate_reference_dict, seq_dict = getAssembleReference(ConnectRefOfEachRefDict)
seq_dict_list.append(seq_dict)
candidate_reference_dict_list.append(candidate_reference_dict)
HGT_ReplicatingRefIndexList, non_HGT_ReplicatingRefIndexList = detectReplicatingRef(candidate_reference_dict, seq_dict, sample_id)
replication_list.append([HGT_ReplicatingRefIndexList, non_HGT_ReplicatingRefIndexList])
genome_sequence_coverage = {}



for i in range(len(seq_dict_list)):
    if args["o"][-1] == '/':
        ref_file = args["o"] + sample_id + '.fasta'
    else:
        ref_file = args["o"] + '/' + sample_id + '.fasta'
    for key in seq_dict_list[i]:
        for j in range(len(seq_dict_list[i][key])):
            ref_name = sample_id + '_' + key + '_hgt_' + str(j)
            ref_seq = seq_dict_list[i][key][j]
            f = open(ref_file, 'a')
            f.write(">" + ref_name + '\n')
            start = 0
            end = 60
            while end < len(ref_seq):
                f.write(ref_seq[start:end]+'\n')
                start = end
                end += 60
            f.write(ref_seq[start:len(ref_seq)]+'\n')
            f.close()

for i in range(len(candidate_reference_dict_list)):
    if args["o"][-1] == '/':
        ref_file = args["o"] + sample_id + '.hgt_reference_info.txt'
    else:
        ref_file = args["o"] + '/' + sample_id + '.hgt_reference_info.txt'
    fo = open(ref_file, "w")
    for key in candidate_reference_dict_list[i]:
        for j in range(len(candidate_reference_dict_list[i][key])):
            ref_name = sample_id + '_' + key + '_hgt_' + str(j)
            reference = candidate_reference_dict_list[i][key][j]
            fo.write(ref_name + '\n')
            for segment in reference:
                seg_info = segment[0] + ',' + str(segment[1]) + ',' + str(segment[2]) + ',' + str(segment[3])
                fo.write(seg_info + '\n')
    fo.close()


for i in range(len(replication_list)):
    if args["o"][-1] == '/':
        output_file = args["o"] + sample_id + '.replication_ref_info.txt'
    else:
        output_file = args["o"] + '/' + sample_id + '.replication_ref_info.txt'
    fo = open(output_file, "w")
    replication_ref_info_list = replication_list[i][0]
    for replication_ref_info in replication_ref_info_list:
        fo.write(sample_id + '_' + replication_ref_info[0] + '_hgt_' + str(replication_ref_info[1]) + ',' + str(replication_ref_info[2]) + '\n')
    replication_ref_info_list = replication_list[i][1]
    for replication_ref_info in replication_ref_info_list:
        fo.write(sample_id + '_' + replication_ref_info[0] + ',' + str(replication_ref_info[1]) + '\n')
    fo.close()
