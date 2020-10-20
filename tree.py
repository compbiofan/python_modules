import copy
from ete3 import Tree
import re
import numpy as np
import sys
sys.path.append('/gpfs/home/xfan2/github/SingleCellCNABenchmark')
import gen_tree
from gen_tree import gen_tree
from cna import check_overlap_one_node, total_edge_length, sep_chromosome_focal_CNA, get_overlap, get_overlap_two_sets
verbose = False

class node():
    def __init__(self, id, parent=-1, children = [], is_leaf = False):
        self.id = id        
        # the edge connecting it with its parent node, containing the bin ids that have the copy number change
        self.edge = []
        self.parent = parent
        self.children = children
        self.is_leaf = is_leaf
        # this entry is nonempty only when it is a leaf
        self.cells = []
        self.level = 0
        # leaf descendant
        self.leaf_desc = []

def get_leaf(tree):
    l = []
    for i in range(len(tree)):
        if tree[i].is_leaf:
            l.append(i)
    return l

def check_circle(tree, node):
    p = tree[node].parent
    ps = [node]
    while p != -1:
        if p in ps:
            print("There is a circle")
            return True
        ps.append(p)
        p = tree[p].parent
    return False 


# given a tree and two nodes, find the lca of the two
def find_lca(tree, node1, node2):
    p = tree[node1].parent
    node1_ps = [node1]
    node2_ps = [node2]
    if check_circle(tree, node1) or check_circle(tree, node2):
        return -1
    while p != -1:
        node1_ps.append(p)
        p = tree[p].parent
    p = node2
    while p not in node1_ps and p != -1:
        node2_ps.append(p)
        p = tree[p].parent
    if p == -1:
        print("Warning: LCA is -1")
        print(node1_ps)
        print(node2_ps)
        print("LCA of " + str(node1) + " and " + str(node2) + " is " + str(p))
    return p 

# check if node1 is node2's ancestor
def if_ancestor(tree, node1, node2):
    p = node2
    while p != -1 and p != node1:
        p = tree[p].parent
    return p == node1

# all events on top of a certain node1 but below node2
def aggregate_events(tree, node1, node2):
    if node1 == node2:
        return []
    e = copy.deepcopy(tree[node1].edge)
    p = tree[node1].parent
    while p != -1 and p != node2:
        if verbose:
            print(str(p) + " edges are ", end = "")
        for i in tree[p].edge:
            e.append(i)
            if verbose:
                print(i, end="")
        if verbose:
            print("")
        p = tree[p].parent
    return e

# given a tree, add the levels to each node, starting from root at level 0
# the tree is from leaves to root
def add_level_reverse(t):
    for i in t[::-1]:
        id = i.id
        if i.parent == -1:
            continue
        else:
            # get the parent's level
            i.level = p_level + 1
    return t 

# this function add the leaf to desc on any tree order (random also permitted)
def add_leaf_desc_general(tree, root):
    if tree[root].is_leaf:
        tree[root].leaf_desc = [root]
        return [root]
    children = tree[root].children
    tree[root].leaf_desc = []
    for i in children:
        leaf_desc = add_leaf_desc_general(tree, i)
        for j in leaf_desc:
            tree[root].leaf_desc.append(j)
    return tree[root].leaf_desc

# for every internal node, add the leaf descendants in leaf_desc
# tree is in reverse order, i.e., leaves' IDs are smaller than internal nodes'
def add_leaf_desc(tree, tree_order):
    rg = range(1, len(tree))
    if tree_order == "forw":
        rg = rg[:,:,-1] 

    for i in rg:
        id = i
        # start from empty before adding
        tree[i].leaf_desc = []
        if tree[i].is_leaf:
            tree[i].leaf_desc.append(id)
        else:
            for j in tree[i].children:
                for x in tree[j].leaf_desc:
                    tree[i].leaf_desc.append(x)

    return tree

# given a binary tree and a node, report the sibling id 
def find_sibling(tree, node):
    cs = tree[tree[node].parent].children
    for i in cs:
        if i != node:
            return i
    return "NA"

# return the robinson-foulds distance between two trees
# the two trees are in newick strings
def rf_dist(t1, t2):
    tree1 = Tree(t1)
    tree2 = Tree(t2)
    r = tree1.compare(tree2)
    return r['norm_rf']

# genenerate newick formatted tree
def gen_newick(Tree, root, option):
    if option == "topdown":
        return convert2newick_topdown(Tree, "(" + str(root) + ")", 0)
    elif option == "bottomup" or option == "general":
        str_ = convert2newick_general(Tree, "(" + str(root) + ")", [root])
        return str_
        

# the leaves in Tree are in the subsequent order as those in t_MyNode
# output a newick string of Tree, with the leaves' id corresponding to t_MyNode
def convert2newick_general(Tree, str_, bucket):
    if len(bucket) == 0:
        return str_
    splitted = re.split(r'(\d+)', str_) 
    new_bucket = []
    for i in bucket:
        for j in range(len(splitted)):
            if str(i) == splitted[j]:
                # replace it with children
                splitted[j] = "(" + ",".join([str(x) for x in Tree[i].children]) + ")"
                for k in Tree[i].children:
                    if not Tree[k].is_leaf:
                        new_bucket.append(k)
    
    return convert2newick_general(Tree, "".join(splitted), new_bucket)
    
# augment newick generation function
def convert2newick_topdown(Tree, str_, i):
    # stop
    if i >= len(Tree):
        return str_
    node = Tree[i]
    parent = node.parent
    children = node.children
    # get rid of the node itself
    c_arr = []
    for c in children:
        if c != str(i):
            c_arr.append(c)
    children = c_arr

    if len(children) >= 1:
        #print children
        # prepare for where to replace
        splitted = re.split(r'(\d+)', str_)
        i_ = -1
        for i_ in range(len(splitted)):
            if str(splitted[i_]) == str(i):
                break
        # now make a new string to replace
        c_ = [str(x) for x in children]
        new_substr = "(" + "," . join(c_) + ")"
        #print str(splitted[i_]) + " is to be replaced by " + new_substr
        # concatenate
        splitted[i_] = new_substr + splitted[i_]
        str_ = "".join(splitted)
        str_ = convert2newick_topdown(Tree, str_, i + 1)
    else:
        str_ = convert2newick_topdown(Tree, str_, i + 1)
    return str_
        
# to be compatible with the old MyNode structure
# add children and identify the is_leaf status
def add_children_from_MyNode(tree_file):
    # add children
    tree = np.load(tree_file, allow_pickle=True)

    # a new tree with a simplified structure
    Tree = []

    # read the npy file and put them in MyNode
    for i in range(len(tree)):
        node_ = tree[i]
        parent = node_.parentID
        ID = i
        Tree.append(node(ID, parent, [], False))
        # update the children of its parent
        if parent != -1:
            Tree[parent].children.append(ID)

    # identify the leaves
    for i in range(len(Tree)):
        if len(Tree[i].children) == 0:
            Tree[i].is_leaf = True 

    return Tree

# from matrix leafid to original leaf id
def reconstruct_leaf_id(leaf_id_f):
    dic_mat2tree = {}
    dic_tree2mat = {}
    f = open(leaf_id_f, "r")
    line = f.readline().rstrip("\n")
    n = 0
    while(line != ""):
        n += 1
        dic_mat2tree[n] = line
        dic_tree2mat[line] = n
        line = f.readline().rstrip("\n")

    return dic_mat2tree, dic_tree2mat 

# replace all ids in newick with their mapped ids
def map_leafid_newick(str_, map_):
    splitted = re.split(r'(\d+)', str_)
    for i in range(len(splitted)):
        if splitted[i].isdigit() and int(splitted[i]) in map_.keys():
            splitted[i] = map_[int(splitted[i])]
    return "".join(splitted)
    
# This function was a modified version of that in read_tree.py in SingleCellBenchmark project. The difference is that instead of output the new CNAs for each node, it creates a new tree (consistent with the tree defined in this file), and add the edges where the copy number chnage occurs. Note that the edge IDs are b followed by the index of the bin, starting from b1. THus the input file segcopy determines the coordinates of the bins.  
# given a segcopy file with all nodes, a tree with the parent-children relationship, output the new CNAs of a child compared to the parent, including internal nodes
def retrieve_new_overlappingCNAs(segcopy_f, Tree, prefix_len):
    # convert it to the tree format in this file (with children and is_leaf and edge)
    my_tree = add_children_from_MyNode(Tree)
    f = open(segcopy_f, "r")
    line = f.readline().rstrip("\n")
    # from leaf id to the column index
    names_h = {}
    # from column index to leaf id
    names_rev = {}
    first = True
    pre_cna = []
    bin_id = 0
    while(line != ""):
        array = re.split(r'\s+', line)
        if first:
            names = array[3:]
            for i in range(len(names)):
                if names[i] == "":
                    break
                name = int(names[i][prefix_len:])
                names_h[name] = i
                names_rev[i] = name
            first = False
        else:
            bin_id += 1
            cnas = array[3:]
            for i in range(len(cnas)):
                if len(pre_cna) == 0:
                    pre_cna = copy.deepcopy(cnas)
                    continue
                if cnas[i] == "":
                    break
                this_id = names_rev[i]
                root = False
                if my_tree[this_id].parent == -1:
                    root = True
                else:
                    p = names_h[my_tree[this_id].parent]
                # check if it is the same as the previous location
                if root and int(cnas[i]) != 2 or int(cnas[i]) - int(pre_cna[i]) != int(cnas[p]) - int(pre_cna[p]):
                    my_tree[this_id].edge.append("b" + str(bin_id))
            pre_cna = copy.deepcopy(cnas)
        line = f.readline().rstrip("\n")

    return my_tree
    
def calculate_overlap(tree, leaves):
    ret_n = []
    ret_l = []
    ret_np = []
    for i in leaves:
        es = tree[i].edge
        n, l, np = get_overlap(es) 
        # overlapping CNA number
        ret_n.append(n)
        # overlapping CNA total length
        ret_l.append(l)
        # percentage of overlapping CNA
        ret_np.append(np)
    return ret_n, ret_l, ret_np

def edge_statistics(tree, leaves):
    ret_t_e = []
    ret_e_avg = []
    for i in leaves:
        n_e = len(tree[i].edge)
        total_l = total_edge_length(tree[i].edge)
        # total CNA number
        ret_t_e.append(n_e)
        # average CNA length
        ret_e_avg.append(float(total_l) / n_e)
    return ret_t_e, ret_e_avg

def edge_statistics_chro_focal(tree, leaves, ref_fa_f, t):
    ret_chro_n = []
    ret_nonchro_n = []
    ret_chro_l_avg = []
    ret_nonchro_l_avg = []
    ret_nonchro_ov_n = []
    ret_nonchro_ov_l = []
    ret_nonchro_ov_p = []

    # overlapping bases among focals, and focal and chromosomals, but not among chromosomals
    ret_ov_l = []

    for i in leaves:
        # statistics on each cna
        chro, nonchro, chro_e, nonchro_e = sep_chromosome_focal_CNA(tree[i].edge, ref_fa_f, t)
        ret_chro_n.append(len(chro))
        ret_nonchro_n.append(len(nonchro))
        ret_chro_l_avg.append(np.average(np.array(chro)))
        ret_nonchro_l_avg.append(np.average(np.array(nonchro)))
    
        # statistics on overlapping can
        #n_chro, l_chro, p_chro = get_overlap(chro_e)
        n_nonchro, l_nonchro, p_nonchro = get_overlap(nonchro_e)
        #ret_chro_ov_n.append(n_chro)
        #ret_chro_ov_l.append(l_chro)
        ret_nonchro_ov_n.append(n_nonchro)
        ret_nonchro_ov_l.append(l_nonchro)
        ret_nonchro_ov_p.append(p_nonchro)

        # for those duplicated chromosomal events, don't count them duplicatedly
        if_merge_chrom = True
        ov_l = get_overlap_two_sets(chro_e, nonchro_e, if_merge_chrom)
        ov_l += l_nonchro
        ret_ov_l.append(ov_l)

    return ret_chro_n, ret_nonchro_n, ret_chro_l_avg, ret_nonchro_l_avg, ret_nonchro_ov_n, ret_nonchro_ov_l, ret_nonchro_ov_p, ret_ov_l



   
# add parent edges to each node that are desc of parent
def add_ancestor_edges(tree, root):
    for i in tree[root].children:
        for j in tree[root].edge:
            tree[i].edge.append(j) 
        add_ancestor_edges(tree, i)
