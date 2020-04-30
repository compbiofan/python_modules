import numpy as np
import sys
import copy
import re

from fa import read_ref_fai_len, de_chr

def total_edge_length(edges):
    total_edge_l = 0
    for i in range(len(edges)):
        cna_l = get_cna_l(edges[i])
        total_edge_l += cna_l
    return total_edge_l

def get_cna_l(e):
    c, s, e = interpret_CNA(e)
    return e - s
        
def check_overlap_one_node(edges):
    n = 0
    l = 0
            
    for i in range(len(edges) - 1):
        for j in range(i+1, len(edges)):
            len_ = if_overlap(edges[i], edges[j])
            if len_ != 0:
                n += 1
                l += len_
    return n, l
        
# chr:s-e
def if_overlap(cna1, cna2):
    c1, s1, e1 = interpret_CNA(cna1)
    c2, s2, e2 = interpret_CNA(cna2)
    l = 0
    if c1 == c2:
        if s1 < e2 and s1 >= s2 or s2 < e1 and s2 > s1:
            if s1 < e2 and s1 >= s2:
                l = min(e2 - s1, e1 - s1)
            else:
                l = min(e1 - s2, e2 - s2)
    return l

# get the overlapping interval 
def overlap_interval(cna1, cna2):
    c1, s1, e1 = interpret_CNA(cna1)
    c2, s2, e2 = interpret_CNA(cna2)
    if s1 < e2 and s1 >= s2:
        return str(s1) + "." + str(min(e2, e1))
    elif s2 < e1 and s2 > s1:
        return str(s2) + "." + str(min(e2, e1))
    else:
        return "NA"
            
# chr:s-e
def interpret_CNA(e):
    a, b = re.split(r':', e)
    c, d = re.split(r'-', b)
    return a, int(c), int(d)
 
# according to the length of the chromosomes shown in ref_fa file, for any CNA > t / chr_len, report it as a chromosomal event. Otherwise it is focal. Return two arrays of the size of the events for each group.  
def sep_chromosome_focal_CNA(edges, ref_fa_f, t):
    c_len = read_ref_fai_len(ref_fa_f)
    chro = []
    nonchro = []
    chro_e = []
    nonchro_e = []
    for i in edges:
        c, s, e = interpret_CNA(i)
        l = e - s
        c = de_chr(c)
        if c in c_len.keys() and l > c_len[c] * t:
            # chromosomal
            chro.append(l)
            chro_e.append(i)
        else:
            nonchro.append(l)
            nonchro_e.append(i)
    return chro, nonchro, chro_e, nonchro_e

# given two sets of cnas, return the overlapping base number between the two sets
def get_overlap_two_sets(chrom_e, focal_e, if_merge_chrom):
    if if_merge_chrom:
        chrom_e = merge_chrom(chrom_e)
        
    overlap_l = 0
    for i in chrom_e:
        for j in focal_e:
            overlap_l += if_overlap(i, j)

    return overlap_l

# merge the same CNAs, return a new array
def merge_chrom(e):
    ret_e = []
    dic_ = {}
    for i in e:
        dic_[i] = 1
    for i in dic_.keys():
        ret_e.append(i)
    return ret_e 
    
# in one array, return the overlapping number, overlapping bases and the percentage of ovlerapping number CNAs 
def get_overlap(es):
    n_e = len(es)
    n_p = n_e * (n_e - 1) / 2
    n, l = check_overlap_one_node(es)
    return n, l, n/float(n_p)


