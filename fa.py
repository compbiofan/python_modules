import numpy as np
import sys
import re

# read a reference fa file, in a dic whose keys are chromosomes (no chr or Chr), store the length of the chromosomes. 
def read_ref_fai_len(fai_f):
    f = open(fai_f, "r")
    l = f.readline()
    ret_dic = {}
    while(l != ""):
        a = re.split(r'\s+', l)
        c_name = a[0]
        c_name = de_chr(c_name)
        c_len = int(a[1])
        ret_dic[c_name] = c_len
        l = f.readline()
    f.close()
    return ret_dic

def de_chr(c):
    if c[0:3] == "chr" or c[0:3] == "CHR" or c[0:3] == "Chr":
        return c[3:]
    else:
        return c
