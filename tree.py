import copy
from ete3 import Tree
import re

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

# all events on top of a certain node
def aggregate_events(tree, node):
    e = copy.deepcopy(tree[node].edge)
    p = tree[node].parent
    while p != -1:
        for i in tree[p].edge:
            e.append(i)
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
        

# the leaves in Tree are in the subsequent order as those in t_MyNode
# output a newick string of Tree, with the leaves' id corresponding to t_MyNode
def convert2newick_general(Tree, str_, bucket):
    if len(bucket) == 0:
        return str_
    splitted = re.split(r'(\d+)', str_) 
    new_bucket = []
    for i in bucket:
        for j in range(len(splitted)):
            if i == int(splitted[j]):
                # replace it with children
                splitted[j] = "(" + ",".join(map(str(x) for x in Tree[i].children)) + ")"
                for k in Tree[i].children:
                    if not Tree[k].is_leaf:
                        new_bucket.append(k)
    
    convert2newick_general(Tree, "".join(splitted), new_bucket)
    
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
        new_substr = "(" + "," . join(children) + ")"
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
    tree = numpy.load(tree_file, allow_pickle=True)

    # a new tree with a simplified structure
    Tree = []

    # read the npy file and put them in MyNode
    for i in range(len(tree)):
        node = tree[i]
        parent = node.parentID
        ID = i
        Tree.append(node(ID, parent, [], False))
        # update the children of its parent
        Tree[parent].children.append(ID)

    # identify the leaves
    for i in range(len(Tree)):
        if len(Tree[i].children) == 0:
            Tree[i].is_leaf = True 

    return Tree

# from matrix leafid to original leaf id
def reconstruct_leaf_id(leaf_id_f):
    dic = []
    f = open(leaf_id_f, "r")
    line = f.readline().rstrip("\n")
    n = 0
    while(line != ""):
        n += 1
        dic[n] = line
        line = f.readline().rstrip("\n")

    return dic 

# replace all ids in newick with their mapped ids
def map_leafid_newick(str_, map_):
    splitted = re.split(r'(\d+)', str_)
    for i in range(len(splitted)):
        if splitted[i] in map_.keys():
            splitted[i] = map_[splitted[i]]
    return "".join(splitted)
    
