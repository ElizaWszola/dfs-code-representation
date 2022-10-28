# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 14:27:16 2021

@author: Chris Wendler
"""
from collections import defaultdict
import numpy as np
from functools import cmp_to_key
from copy import deepcopy

# https://sites.cs.ucsb.edu/~xyan/papers/gSpan.pdf
def edgeLt(e1, e2):
    if e1[0] < e1[1] and e2[0] < e2[1]:
        return e1[1] < e2[1]
    elif e1[0] > e1[1] and e2[0] > e2[1]:
        return (e1[0] < e2[0]) or (e1[0] == e2[0] and e1[1] < e2[1])
    elif e1[0] < e1[1] and e2[0] > e2[1]:
        return e1[1] <= e2[0]
    else:
        return e1[0] < e2[1]

# https://sites.cs.ucsb.edu/~xyan/papers/gSpan.pdf there is the caveat: E_T1 yields a total order on the edges
# E_T2 yields another total order on the edges... therefore we just take this definition from the paper
def tupleLt(a, b):
    if a[0] > a[1] and b[0] < b[1]:
        return True
    elif a[0] > a[1] and b[0] > b[1] and a[1] < b[1]:
        return True
    elif a[0] > a[1] and b[0] > b[1] and a[1] == b[1] and ((a[2], a[3], a[4]) < (b[2], b[3], b[4])):
        return True
    elif a[0] < a[1] and b[0] < b[1] and a[0] > b[0]:
        return True
    elif a[0] < a[1] and b[0] < b[1] and a[0] == b[0] and ((a[2], a[3], a[4]) < (b[2], b[3], b[4])):
        return True
    else:
        return False


def codeLt(c1, c2):
    m = len(c1)
    n = len(c2)
    for idx, (t1, t2) in enumerate(zip(c1[:min(m, n)], c2[:min(m, n)])):
        #print(t1 != t2, tupleLt(t1, t2), t1, t2)
        if tuple(t1[:-3]) != tuple(t2[:-3]):
            return tupleLt(t1, t2)
    return m < n


def cmp(e1, e2):
    if e1 == e2:
        return 0
    elif edgeLt(e1, e2):
        return -1
    else:
        return 1


def isValidDFSCode(code):
    res = True
    for i in range(1, len(code)):
        e1 = code[i-1][:2]
        e2 = code[i][:2]
        res = res and edgeLt(e1, e2)
    return res     


class Graph:
    """
    Starting point: https://www.geeksforgeeks.org/depth-first-search-or-dfs-for-a-graph/
    """
 
    def __init__(self, edge_list, vertex_labels, edge_labels):
        self.graph = defaultdict(list)
        self.edge_list = edge_list
        self.vertex_labels = vertex_labels
        self.edge_labels = edge_labels
        self.edge_index = {}
        self.edge_index_undirected = {}
        for idx, edge in enumerate(edge_list):
            self.addEdge(edge[0], edge[1])
            self.edge_index[(edge[0], edge[1])] = idx
            if (edge[0], edge[1]) not in self.edge_index_undirected:
                eidx = len(self.edge_index_undirected)//2
                self.edge_index_undirected[(edge[0], edge[1])] = eidx
                self.edge_index_undirected[(edge[1], edge[0])] = eidx
 
    
    def addEdge(self, u, v):
        self.graph[u].append(v)
        
        
    def DFSUtil(self, v, visited, dfs_index):
        # Mark the current node as visited
        # and print it
        dfs_index[v] = len(visited)
        visited.add(v)        
 
        # Recur for all the vertices
        # adjacent to this vertex
        for neighbour in np.random.permutation(self.graph[v]):
            if neighbour not in visited:
                self.DFSUtil(neighbour, visited, dfs_index)
            
 
    def DFS(self, v):
        # Create a set to store visited vertices
        visited = set()
        dfs_index = {}
        # Call the recursive helper function
        # to print DFS traversal
        self.DFSUtil(v, visited, dfs_index)
        return dfs_index

    
    def DFSCodeUtil(self, u, v, evisited, dfs_index, code):
        vlabel = self.vertex_labels
        elabel = self.edge_labels 
        dfs_index[v] = len(dfs_index)
        eidx = self.edge_index[(u, v)]
        code += [(dfs_index[u], dfs_index[v], vlabel[u], elabel[eidx], vlabel[v], u, eidx, v)]
        #print(code)
        #print("+++++++++++++++")
        evisited.add(self.edge_index_undirected[(u, v)])
        
        # add backward edges
        bwd_targets = {dfs_index[neighbour] : neighbour for neighbour in self.graph[v] if neighbour in dfs_index}
        for dfs_to in sorted(bwd_targets.keys()):
            u = bwd_targets[dfs_to]
            eidx = self.edge_index[(v, u)]
            if self.edge_index_undirected[(v, u)] not in evisited:
                code += [(dfs_index[v], dfs_to, vlabel[v], elabel[eidx], vlabel[u], v, eidx, u)]
                evisited.add(self.edge_index_undirected[(v, u)])
        # add forward edges
        for neighbour in np.random.permutation(self.graph[v]):
            if neighbour not in dfs_index:
                self.DFSCodeUtil(v, neighbour, evisited, dfs_index, code) 
                
    
    def DFSCode(self):
        """

        Returns
        -------
        code : dfs code, where each entry is made up from (dfs_from, 
                                                           dfs_to,
                                                           feat_from,
                                                           feat_edge,
                                                           feat_to,
                                                           idx_from,
                                                           idx_edge,
                                                           idx_to),
        here idx_edge is the index of the edge in the edge list (in contrast,
        the cpp implementation of minimal dfs codes puts the undirected edge index here).
        This is taken care of in utils.py.
        dfs_index: list where entry i corresponds to dfs_index[i]
        """
        #print(self.edge_list)
        start = np.random.randint(len(self.vertex_labels))
        dfs_index = {start: 0}
        evisited = set()
        code = []
        for to in np.random.permutation(self.graph[start]):
            if to not in dfs_index:
                self.DFSCodeUtil(start, to, evisited, dfs_index, code) 
        if len(code) == 0: # deal with single atom molecules...
            self.DFSCodeUtil(start, self.graph[start][0], evisited, dfs_index, code) 
        return code, [dfs_index[u] for u in range(len(self.vertex_labels))]
    
              
    def MinDFSCodeUtil(self, u, v, evisited, dfs_index, code, rightmost_path, codeLt = codeLt):
        # update data structures by the new edge
        vlabel = self.vertex_labels
        elabel = self.edge_labels 
        dfs_index[v] = len(dfs_index)
        eidx = self.edge_index[(u, v)]
        evisited.add(self.edge_index_undirected[(u, v)])
        rightmost_path += [v]
        code += [(dfs_index[u], dfs_index[v], vlabel[u], elabel[eidx], vlabel[v], u, self.edge_index_undirected[(u, v)], v)]
        #print(code)
        if not codeLt(code, self.min_dfs_code) and len(self.min_dfs_code) != 0:
            return
        # because we consider undirected graphs all backward edges must target a
        # vertex on the rightmost path
        bwd_targets = {dfs_index[w] : w for w in self.graph[v] if w in rightmost_path}
        for dfs_to in sorted(bwd_targets.keys()):
            w = bwd_targets[dfs_to]
            eidx = self.edge_index[(v, w)]
            if self.edge_index_undirected[(v, w)] not in evisited:
                code += [(dfs_index[v], dfs_to, vlabel[v], elabel[eidx], vlabel[w], v, self.edge_index_undirected[(v, w)], w)]
                evisited.add(self.edge_index_undirected[(v, w)])
                if not codeLt(code, self.min_dfs_code) and len(self.min_dfs_code) != 0:
                    return
        
        # if we are done we overwrite
        if len(code) == len(self.edge_list)//2:
            self.min_dfs_code = code
            self.min_dfs_index = dfs_index
            #print("_________________________________________________")
                
        # add forward edges in a way that is consistent with our DFS implementation
        # first we try to add the ones at the rightmost vertex v
        # second we try to add the ones on the rightmost path (starting from 
        # vertices that are close to v going all the way up to the root)
        for idx, v in enumerate(rightmost_path[::-1]):
            ranks = np.unique([self.edge_rank[(vlabel[v], 
                                               elabel[self.edge_index[(v, w)]], 
                                               vlabel[w])] 
                               for w in self.graph[v]
                               if self.edge_index_undirected[(v, w)] not in evisited and w not in dfs_index])
            for w in self.graph[v]:
                eidx = self.edge_index[(v, w)]
                label = (vlabel[v], elabel[eidx], vlabel[w])
                if self.edge_index_undirected[(v, w)] not in evisited and \
                   (w not in dfs_index) and (self.edge_rank[label] == ranks[0]):
                    self.MinDFSCodeUtil(v, w, deepcopy(evisited),
                                        deepcopy(dfs_index), 
                                        deepcopy(code),
                                        deepcopy(rightmost_path[:len(rightmost_path)-idx])) 
        
        
    def MinDFSCode(self, codeLt = codeLt):
        #print("------------------------------------------------------")
        vlabel = self.vertex_labels
        elabel = self.edge_labels
        # 1. rank edges
        unique_edge_labels = sorted(set([(vlabel[e[0]], elabel[self.edge_index[(e[0], e[1])]], vlabel[e[1]]) for e in self.edge_list]))
        self.edge_rank = {label:rank for rank, label in enumerate(unique_edge_labels)}
        self.min_dfs_code = []
        self.min_dfs_index = {}
        # 2. start one computation from each minimal edge
        for e in self.edge_list:
            label = (vlabel[e[0]], elabel[self.edge_index[(e[0], e[1])]], vlabel[e[1]])
            if self.edge_rank[label] == 0:
                self.MinDFSCodeUtil(e[0], e[1], set(), {e[0]: 0}, [], [e[0]], codeLt = codeLt)
        return self.min_dfs_code, [self.min_dfs_index[u] for u in range(len(self.vertex_labels))]
        
        
        
    def DFSCodeLegacy(self):
        start = np.random.randint(len(self.vertex_labels))
        dfs_index = self.DFS(start)
        dfs2edge = {(dfs_index[e[0]], dfs_index[e[1]]):e for e in self.edge_list}
        dfs_edges = sorted(list(dfs2edge.keys()), key=cmp_to_key(cmp))
        edge_history = []
        while len(dfs_edges) > 0:
            for e in dfs_edges:
                if len(edge_history) == 0 or edgeLt(edge_history[-1], e):
                    break
            edge_history += [e]
            dfs_edges.remove(e)
            try:
                dfs_edges.remove((e[1], e[0]))
            except ValueError:
                continue
        edge_history = [dfs2edge[e] for e in edge_history]
        
        code = []
        vlabel = self.vertex_labels
        elabel = self.edge_labels 
        for (u, v) in edge_history:
            eidx = self.edge_index[(u, v)]
            code += [(dfs_index[u], dfs_index[v], vlabel[u], elabel[eidx], vlabel[v], u, eidx, v)]
        return code, [dfs_index[u] for u in range(len(self.vertex_labels))]
        