#!/home/xf2/.conda/envs/my_root/bin/python
import sys
import numpy as np
from man_file import read_csv_f, print_matrix

if len(sys.argv) <= 1:
    print("""
    Given a folder with wig.csv files, a leaf node ID (one per row) list file, output a SegCopy formatted file that has the cells as the columns, the rows the genomic regions, the entries the read counts. 
    Usage: python
    """ + sys.argv[0] + """ [wigcsv_folder] [leaf.id] [prefix_before_id] [postfix_after_id] [output_segcopy_f] [printing_option:1(without row and column header);2(with both)""")
    sys.exit(0)

folder = sys.argv[1]
leaf_f = sys.argv[2]
prefix = sys.argv[3]
postfix = sys.argv[4]
output_f = sys.argv[5]
option = sys.argv[6]
if option == '2':
    option = True
else:
    option = False

# open the leaf list file
sep = ","
header = True
f = open(leaf_f, "r")
line = f.readline().rstrip("\n")
m = []
genomic_a = []
leaf_a = []
while(line != ""):
    leaf_id = line
    csv_f = folder + "/" + prefix + leaf_id + postfix
    # csv file is separated 
    a = read_csv_f(csv_f, sep, header, [6])
    m.append(a)
    leaf_a.append(line)
    if len(genomic_a) == 0:
        genomic_a = read_csv_f(csv_f, sep, header, [0, 1, 2])
    line = f.readline().rstrip("\n")

mT = np.array(m).T
    
# column, row
print_matrix(mT, option, genomic_a, leaf_a, output_f)
