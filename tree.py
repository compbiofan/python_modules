import copy
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
