#!/usr/bin/python

# module for union find to get connected component

# get the root given a node and a tree
def root(p, tr):
    while tr[p] != p:
        p = tr[p]
    return p

# connect the two trees, smaller to the larger
def connect(p, q, tr, sz):
    if sz[p] > sz[q]:
        tr[q] = p
        sz[p] += sz[q]
    else:
        tr[p] = q
        sz[q] += sz[p]
    return (tr, sz)

