import unittest
import numpy as np
import _dfs_codes as dfs
import dfs_code

import networkx as nx
import tqdm
import time
from scipy.io import loadmat

n_tries = 100
n_vertices = 10
p_edge = 0.3
max_labels = 15
timeout = int(120*1e6)

def adjacency2list(A, node_labels):
    edge_list = np.asarray(np.where(A)).T.tolist()
    idx_dict = {}
    idx = 0
    for edge in edge_list:
        key = tuple([edge[1], edge[0]])
        edge_label = node_labels[edge[0]]*node_labels[edge[1]]
        if key in idx_dict:
            edge += [edge_label, idx_dict[key]]
        else:
            idx_dict[tuple(edge)] = idx
            idx += 1
            edge += [edge_label, idx_dict[tuple(edge)]]
    return edge_list

# https://sites.cs.ucsb.edu/~xyan/papers/gSpan.pdf
def edgeLeq(e1, e2):
    if e1[0] > e1[1] and e2[0] < e2[1] and e1[0] < e2[1]:
        return True
    if e1[0] < e1[1] and e2[0] > e2[1] and e1[1] <= e2[0]:
        return True
    if e1[0] < e1[1] and e2[0] < e2[1] and e1[1] < e2[1]:
        return True
    if e1[0] > e1[1] and e2[0] > e2[1] and (e1[0] < e2[0] or (e1[0] == e2[0] and e1[1] < e2[1])):
        return True
    return False    
    # this is garbage from : https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1184038
    #if e1[0] == e2[0] and e1[1] < e2[1]:
    #    return True
    #elif e1[0] < e1[1] and e1[1] == e2[0]:
    #    return True
    #return False

def isValidDFSCode(code):
    res = True
    for i in range(1, len(code)):
        e1 = code[i-1][:2]
        e2 = code[i][:2]
        res = res and edgeLeq(e1, e2)
        
    return res        

class MinDFSTest(unittest.TestCase):   
    def test_random_dfs_codes(self):
        print('testing random dfs codes...')
        for tries in tqdm.tqdm(range(n_tries)):
            A = np.random.binomial(1, p_edge, (n_vertices, n_vertices))
            A = np.minimum(A + A.T, 1)
            for i in range(A.shape[0]):
                A[i, i] = 0
            n_types = np.random.randint(1, max_labels)
            if n_types > 1:
                node_types = np.random.randint(1, n_types, A.shape[0]).tolist()
            else:
                node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
            edge_list = adjacency2list(A, node_types)
            
            
            g1 = nx.from_numpy_matrix(A)
            
            if not nx.is_connected(g1):
                continue
            
            edge_list = np.asarray(edge_list)
            edges = edge_list[:, :2].tolist()
            elabels = edge_list[:, 2]
            G = dfs_code.Graph(edges, node_types, elabels)
            code, dfs_index = G.DFSCode()
            

            self.assertTrue(len(code) == len(edge_list)//2)
            self.assertTrue(isValidDFSCode(code))
            
            edge_set = []
            for r in code:
                edge_set += [(r[5], r[7]),(r[7], r[5])]
            edge_set = set(edge_set)
                     
		    # check all vertex and edge ids
            for e in edge_list:
                self.assertTrue((e[0], e[1]) in edge_set)

"""
    def test_exchanging_order_on_python_codes(self):
        print('testing python minimum dfs codes with the std. lexicographic order on 5-tuples...')
        def codeLt(c1, c2):
            m = len(c1)
            n = len(c2)
            for idx, (t1, t2) in enumerate(zip(c1[:min(m, n)], c2[:min(m, n)])):
                #print(t1 != t2, tupleLt(t1, t2), t1, t2)
                a = tuple(t1[:-3])
                b = tuple(t2[:-3])
                if a != b:
                    return a < b
            return m < n
            
            
        for tries in tqdm.tqdm(range(n_tries)):
            A = np.random.binomial(1, p_edge, (n_vertices, n_vertices))
            A = np.minimum(A + A.T, 1)
            for i in range(A.shape[0]):
                A[i, i] = 0
            n_types = np.random.randint(1, max_labels)
            if n_types > 1:
                node_types = np.random.randint(1, n_types, A.shape[0]).tolist()
            else:
                node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
            edge_list = adjacency2list(A, node_types)
            
            # permute the vertices and edges accordingly
            perm = np.random.permutation(A.shape[0])
            A_ = A.copy()[perm]
            A_ = A_[:, perm]
            
            g1 = nx.from_numpy_matrix(A)
            g2 = nx.from_numpy_matrix(A_)
            
            if not nx.is_connected(g1):
                continue
            
            self.assertTrue(nx.is_isomorphic(g1, g2))
            node_types_ = np.asarray(node_types).copy()[perm].tolist()
            edge_list_ = adjacency2list(A_, node_types_)
            
            
            edge_list = np.asarray(edge_list)
            edges = edge_list[:, :2].tolist()
            elabels = edge_list[:, 2]
            G = dfs_code.Graph(edges, node_types, elabels)
            code, dfs_index = G.MinDFSCode(codeLt=codeLt)
            rnd, _ = G.DFSCode()
            
            edge_list_ = np.asarray(edge_list_)
            edges_ = edge_list_[:, :2].tolist()
            elabels_ = edge_list_[:, 2]
            G_ = dfs_code.Graph(edges_, node_types_, elabels_)
            code_, dfs_index_ = G_.MinDFSCode(codeLt=codeLt)
            rnd_, _ = G_.DFSCode()
            

            self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
            self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))
            self.assertTrue(len(code) == len(edge_list)//2)
            self.assertTrue(len(code_) == len(edge_list_)//2)
            self.assertTrue(isValidDFSCode(code))
            self.assertTrue(isValidDFSCode(code_))

            
            edge_dict = {r[6]: (r[5], r[7], r[3]) for r in code}  
            edge_dict_ = {r[6]: (r[5], r[7], r[3]) for r in code_}  
            #print(A)
            #print(A_)
            for row, row_ in zip(code, code_):
			    # check dfs code
                #print(row)
                #print(row_)
                #print("--------------------------")
                for el, el_ in zip(row[:5], row_[:5]):
                    self.assertEqual(el, el_)
                
		    # check all vertex and edge ids
            for e, e_ in zip(edge_list, edge_list_):
                eid = e[-1]
                eid_ = e_[-1]
                de = edge_dict[eid]
                self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
                de = edge_dict_[eid_]
                e = e_
                self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
    
    def test_robustness_to_graph_isomorphy_python(self):
        print('testing python minimum dfs codes...')
        for tries in tqdm.tqdm(range(n_tries)):
            A = np.random.binomial(1, p_edge, (n_vertices, n_vertices))
            A = np.minimum(A + A.T, 1)
            for i in range(A.shape[0]):
                A[i, i] = 0
            n_types = np.random.randint(1, max_labels)
            if n_types > 1:
                node_types = np.random.randint(1, n_types, A.shape[0]).tolist()
            else:
                node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
            edge_list = adjacency2list(A, node_types)
            
            # permute the vertices and edges accordingly
            perm = np.random.permutation(A.shape[0])
            A_ = A.copy()[perm]
            A_ = A_[:, perm]
            
            g1 = nx.from_numpy_matrix(A)
            g2 = nx.from_numpy_matrix(A_)
            
            if not nx.is_connected(g1):
                continue
            
            self.assertTrue(nx.is_isomorphic(g1, g2))
            node_types_ = np.asarray(node_types).copy()[perm].tolist()
            edge_list_ = adjacency2list(A_, node_types_)
            
            code_cpp = dfs.compute_minimal_dfs_code(edge_list, node_types, timeout)
            code_cpp_ = dfs.compute_minimal_dfs_code(edge_list_, node_types_, timeout)
            
            edge_list = np.asarray(edge_list)
            edges = edge_list[:, :2].tolist()
            elabels = edge_list[:, 2]
            G = dfs_code.Graph(edges, node_types, elabels)
            code, dfs_index = G.MinDFSCode()
            rnd, _ = G.DFSCode()
            
            edge_list_ = np.asarray(edge_list_)
            edges_ = edge_list_[:, :2].tolist()
            elabels_ = edge_list_[:, 2]
            G_ = dfs_code.Graph(edges_, node_types_, elabels_)
            code_, dfs_index_ = G_.MinDFSCode()
            rnd_, _ = G_.DFSCode()
            
            
            
            self.assertTrue(np.sum(np.abs(np.asarray(code)[:, :5] - np.asarray(code_cpp)[:, :5])) == 0)
            self.assertTrue(np.sum(np.abs(np.asarray(code_)[:, :5] - np.asarray(code_cpp_)[:, :5])) == 0)
            self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
            self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))
            self.assertTrue(len(code) == len(edge_list)//2)
            self.assertTrue(len(code_) == len(edge_list_)//2)
            self.assertTrue(isValidDFSCode(code))
            self.assertTrue(isValidDFSCode(code_))

            
            edge_dict = {r[6]: (r[5], r[7], r[3]) for r in code}  
            edge_dict_ = {r[6]: (r[5], r[7], r[3]) for r in code_}  
            #print(A)
            #print(A_)
            for row, row_ in zip(code, code_):
			    # check dfs code
                #print(row)
                #print(row_)
                #print("--------------------------")
                for el, el_ in zip(row[:5], row_[:5]):
                    self.assertEqual(el, el_)
                
		    # check all vertex and edge ids
            for e, e_ in zip(edge_list, edge_list_):
                eid = e[-1]
                eid_ = e_[-1]
                de = edge_dict[eid]
                self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
                de = edge_dict_[eid_]
                e = e_
                self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
    
    def test_robustness_to_graph_isomorphy(self):
        print('testing minimum dfs codes...')
        for tries in tqdm.tqdm(range(n_tries)):
            A = np.random.binomial(1, p_edge, (n_vertices, n_vertices))
            A = np.minimum(A + A.T, 1)
            for i in range(A.shape[0]):
                A[i, i] = 0
            n_types = np.random.randint(1, max_labels)
            if n_types > 1:
                node_types = np.random.randint(1, n_types, A.shape[0]).tolist()
            else:
                node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
            edge_list = adjacency2list(A, node_types)
            
            # permute the vertices and edges accordingly
            perm = np.random.permutation(A.shape[0])
            A_ = A.copy()[perm]
            A_ = A_[:, perm]
            
            g1 = nx.from_numpy_matrix(A)
            g2 = nx.from_numpy_matrix(A_)
            
            if not nx.is_connected(g1):
                continue
            
            self.assertTrue(nx.is_isomorphic(g1, g2))
            
            node_types_ = np.asarray(node_types).copy()[perm].tolist()
            edge_list_ = adjacency2list(A_, node_types_)
            code = dfs.compute_minimal_dfs_code(edge_list, node_types, timeout)
            code_ = dfs.compute_minimal_dfs_code(edge_list_, node_types_, timeout)
            
            edge_list = np.asarray(edge_list)
            edges = edge_list[:, :2].tolist()
            elabels = edge_list[:, 2]
            G = dfs_code.Graph(edges, node_types, elabels)
            rnd, _ = G.DFSCode()
            
            edge_list_ = np.asarray(edge_list_)
            edges_ = edge_list_[:, :2].tolist()
            elabels_ = edge_list_[:, 2]
            G_ = dfs_code.Graph(edges_, node_types_, elabels_)
            rnd_, _ = G_.DFSCode()
            
            
            self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
            self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))

            self.assertTrue(len(code) == len(edge_list)//2)
            self.assertTrue(len(code_) == len(edge_list_)//2)
            self.assertTrue(isValidDFSCode(code))
            self.assertTrue(isValidDFSCode(code_))
            
            edge_dict = {r[6]: (r[5], r[7], r[3]) for r in code}  
            edge_dict_ = {r[6]: (r[5], r[7], r[3]) for r in code_}  
            for row, row_ in zip(code, code_):
			    # check dfs code
                for el, el_ in zip(row[:5], row_[:5]):
                    self.assertEqual(el, el_)
		    # check all vertex and edge ids
            for e, e_ in zip(edge_list, edge_list_):
                eid = e[-1]
                eid_ = e_[-1]
                de = edge_dict[eid]
                self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
                de = edge_dict_[eid_]
                e = e_
                self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
    
    
    def test_robustness_to_graph_isomorphy_moresophisticated(self):
        print('testing minimum dfs codes...')
        graphs = loadmat("random_graphs/RandomGraphs.mat")
        types = ['barabasiAlbertGraphs', 'klemmEguilezGraphs', 'wattStrogatzGraphs']
        for t in types:
            print('on %s...'%t)
            for A in tqdm.tqdm(graphs[t].T):
                A = np.minimum(A + A.T, 1)
                for i in range(A.shape[0]):
                    A[i, i] = 0
                n_types = np.random.randint(1, max_labels)
                if n_types > 1:
                    node_types = np.random.randint(1, n_types, A.shape[0]).tolist()
                else:
                    node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
                edge_list = adjacency2list(A, node_types)
                
                # permute the vertices and edges accordingly
                perm = np.random.permutation(A.shape[0])
                A_ = A.copy()[perm]
                A_ = A_[:, perm]
                
                g1 = nx.from_numpy_matrix(A)
                g2 = nx.from_numpy_matrix(A_)
                
                if not nx.is_connected(g1):
                    continue
                
                self.assertTrue(nx.is_isomorphic(g1, g2))
                
                node_types_ = np.asarray(node_types).copy()[perm].tolist()
                edge_list_ = adjacency2list(A_, node_types_)
                
                code = dfs.compute_minimal_dfs_code(edge_list, node_types, timeout)
                code_ = dfs.compute_minimal_dfs_code(edge_list_, node_types_, timeout)
                
                edge_list = np.asarray(edge_list)
                edges = edge_list[:, :2].tolist()
                elabels = edge_list[:, 2]
                G = dfs_code.Graph(edges, node_types, elabels)
                rnd, _ = G.DFSCode()
                
                edge_list_ = np.asarray(edge_list_)
                edges_ = edge_list_[:, :2].tolist()
                elabels_ = edge_list_[:, 2]
                G_ = dfs_code.Graph(edges_, node_types_, elabels_)
                rnd_, _ = G_.DFSCode()
                
                self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
                self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))
    
                self.assertTrue(len(code) == len(edge_list)//2)
                self.assertTrue(len(code_) == len(edge_list_)//2)
                self.assertTrue(isValidDFSCode(code))
                self.assertTrue(isValidDFSCode(code_))
                
                edge_dict = {r[6]: (r[5], r[7], r[3]) for r in code}  
                edge_dict_ = {r[6]: (r[5], r[7], r[3]) for r in code_}  
                for row, row_ in zip(code, code_):
    			    # check dfs code
                    for el, el_ in zip(row[:5], row_[:5]):
                        self.assertEqual(el, el_)
    		    # check all vertex and edge ids
                for e, e_ in zip(edge_list, edge_list_):
                    eid = e[-1]
                    eid_ = e_[-1]
                    de = edge_dict[eid]
                    self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
                    de = edge_dict_[eid_]
                    e = e_
                    self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
    
    
    def test_robustness_to_graph_isomorphy_python_moresophisticated(self):
        #return
        print('testing python minimum dfs codes...')
        graphs = loadmat("random_graphs/RandomGraphs.mat")
        types = ['barabasiAlbertGraphs', 'klemmEguilezGraphs', 'wattStrogatzGraphs']
        for t in types:
            print('on %s...'%t)
            for A in tqdm.tqdm(graphs[t].T):
                A = np.minimum(A + A.T, 1)
                for i in range(A.shape[0]):
                    A[i, i] = 0
                n_types = np.random.randint(1, max_labels)
                if n_types > 1:
                    node_types = np.random.randint(1, n_types, A.shape[0]).tolist()
                else:
                    node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
                edge_list = adjacency2list(A, node_types)
                
                # permute the vertices and edges accordingly
                perm = np.random.permutation(A.shape[0])
                A_ = A.copy()[perm]
                A_ = A_[:, perm]
                
                node_types_ = np.asarray(node_types)[perm].tolist()
                edge_list_ = adjacency2list(A_, node_types_)
                
                code_cpp = dfs.compute_minimal_dfs_code(edge_list, node_types, timeout)
                code_cpp_ = dfs.compute_minimal_dfs_code(edge_list_, node_types_, timeout)
                
                g1 = nx.from_numpy_matrix(A)
                g2 = nx.from_numpy_matrix(A_)
                
                if not nx.is_connected(g1):
                    continue
                
                self.assertTrue(nx.is_isomorphic(g1, g2))
                
                edge_list = np.asarray(edge_list)
                edges = edge_list[:, :2].tolist()
                elabels = edge_list[:, 2]
                G = dfs_code.Graph(edges, node_types, elabels)
                code, dfs_index = G.MinDFSCode()
                rnd, _ = G.DFSCode()
                
                edge_list_ = np.asarray(edge_list_)
                edges_ = edge_list_[:, :2].tolist()
                elabels_ = edge_list_[:, 2]
                G_ = dfs_code.Graph(edges_, node_types_, elabels_)
                code_, dfs_index_ = G_.MinDFSCode()
                rnd_, _ = G_.DFSCode()
                
                
                self.assertTrue(np.sum(np.abs(np.asarray(code)[:, :5] - np.asarray(code_cpp)[:, :5])) == 0)
                self.assertTrue(np.sum(np.abs(np.asarray(code_)[:, :5] - np.asarray(code_cpp_)[:, :5])) == 0)
                
                self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
                self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))
                
                
                self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
                self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))
    
                self.assertTrue(len(code) == len(edge_list)//2)
                self.assertTrue(len(code_) == len(edge_list_)//2)
                self.assertTrue(isValidDFSCode(code))
                self.assertTrue(isValidDFSCode(code_))
                
                edge_dict = {r[6]: (r[5], r[7], r[3]) for r in code}  
                edge_dict_ = {r[6]: (r[5], r[7], r[3]) for r in code_}  
                for row, row_ in zip(code, code_):
    			    # check dfs code
                    for el, el_ in zip(row[:5], row_[:5]):
                        self.assertEqual(el, el_)
    		    # check all vertex and edge ids
                for e, e_ in zip(edge_list, edge_list_):
                    eid = e[-1]
                    eid_ = e_[-1]
                    de = edge_dict[eid]
                    self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
                    de = edge_dict_[eid_]
                    e = e_
                    self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
                    
                    
    def test_robustness_to_graph_isomorphy_python_moresophisticated(self):
        print('testing python minimum dfs codes with different codeLt...')
        def codeLt(c1, c2):
            m = len(c1)
            n = len(c2)
            for idx, (t1, t2) in enumerate(zip(c1[:min(m, n)], c2[:min(m, n)])):
                #print(t1 != t2, tupleLt(t1, t2), t1, t2)
                a = tuple(t1[:-3])
                b = tuple(t2[:-3])
                if a != b:
                    return a < b
            return m < n
        graphs = loadmat("random_graphs/RandomGraphs.mat")
        types = ['barabasiAlbertGraphs', 'klemmEguilezGraphs', 'wattStrogatzGraphs']
        for t in types:
            print('on %s...'%t)
            for A in tqdm.tqdm(graphs[t].T):
                A = np.minimum(A + A.T, 1)
                for i in range(A.shape[0]):
                    A[i, i] = 0
                n_types = np.random.randint(1, max_labels)
                if n_types > 1:
                    node_types = np.random.randint(1, n_types, A.shape[0]).tolist()
                else:
                    node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
                edge_list = adjacency2list(A, node_types)
                
                # permute the vertices and edges accordingly
                perm = np.random.permutation(A.shape[0])
                A_ = A.copy()[perm]
                A_ = A_[:, perm]
                
                node_types_ = np.asarray(node_types)[perm].tolist()
                edge_list_ = adjacency2list(A_, node_types_)
                
                
                g1 = nx.from_numpy_matrix(A)
                g2 = nx.from_numpy_matrix(A_)
                
                if not nx.is_connected(g1):
                    continue
                
                self.assertTrue(nx.is_isomorphic(g1, g2))
                
                edge_list = np.asarray(edge_list)
                edges = edge_list[:, :2].tolist()
                elabels = edge_list[:, 2]
                G = dfs_code.Graph(edges, node_types, elabels)
                code, dfs_index = G.MinDFSCode(codeLt=codeLt)
                rnd, _ = G.DFSCode()
                
                edge_list_ = np.asarray(edge_list_)
                edges_ = edge_list_[:, :2].tolist()
                elabels_ = edge_list_[:, 2]
                G_ = dfs_code.Graph(edges_, node_types_, elabels_)
                code_, dfs_index_ = G_.MinDFSCode(codeLt=codeLt)
                rnd_, _ = G_.DFSCode()
                
                                
                self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
                self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))
                
                
                self.assertTrue((rnd == code) or dfs_code.codeLt(code, rnd))
                self.assertTrue((rnd_ == code_) or dfs_code.codeLt(code_, rnd_))
    
                self.assertTrue(len(code) == len(edge_list)//2)
                self.assertTrue(len(code_) == len(edge_list_)//2)
                self.assertTrue(isValidDFSCode(code))
                self.assertTrue(isValidDFSCode(code_))
                
                edge_dict = {r[6]: (r[5], r[7], r[3]) for r in code}  
                edge_dict_ = {r[6]: (r[5], r[7], r[3]) for r in code_}  
                for row, row_ in zip(code, code_):
    			    # check dfs code
                    for el, el_ in zip(row[:5], row_[:5]):
                        self.assertEqual(el, el_)
    		    # check all vertex and edge ids
                for e, e_ in zip(edge_list, edge_list_):
                    eid = e[-1]
                    eid_ = e_[-1]
                    de = edge_dict[eid]
                    self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
                    de = edge_dict_[eid_]
                    e = e_
                    self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))

    def test_robustness_to_graph_isomorphy_fixed(self):
        print('testing minimum dfs code for fixed graph...')
        edge_list = [[0, 3, 24, 1], [0, 4, 20, 2], [0, 5, 8, 3], [0, 6, 24, 4], [0, 8, 8, 5], [1, 4, 25, 6], [1, 7, 5, 7], [2, 4, 30, 8], [2, 5, 12, 9], [2, 8, 12, 10], [2, 9, 18, 11], [3, 0, 24, 1], [3, 4, 30, 12], [3, 5, 12, 13], [4, 0, 20, 2], [4, 1, 25, 6], [4, 2, 30, 8], [4, 3, 30, 12], [4, 7, 5, 14], [5, 0, 8, 3], [5, 2, 12, 9], [5, 3, 12, 13], [5, 8, 4, 15], [6, 0, 24, 4], [6, 8, 12, 16], [6, 9, 18, 17], [7, 1, 5, 7], [7, 4, 5, 14], [7, 9, 3, 18], [8, 0, 8, 5], [8, 2, 12, 10], [8, 5, 4, 15], [8, 6, 12, 16], [9, 2, 18, 11], [9, 6, 18, 17], [9, 7, 3, 18]]
        edge_list_ = [[0, 1, 12, 1], [0, 4, 4, 2], [0, 5, 8, 3], [0, 6, 12, 4], [1, 0, 12, 1], [1, 4, 12, 5], [1, 7, 18, 6], [1, 8, 30, 7], [2, 3, 5, 8], [2, 8, 25, 9], [3, 2, 5, 8], [3, 7, 3, 10], [3, 8, 5, 11], [4, 0, 4, 2], [4, 1, 12, 5], [4, 5, 8, 12], [4, 9, 12, 13], [5, 0, 8, 3], [5, 4, 8, 12], [5, 6, 24, 14], [5, 8, 20, 15], [5, 9, 24, 16], [6, 0, 12, 4], [6, 5, 24, 14], [6, 8, 30, 17], [7, 1, 18, 6], [7, 3, 3, 10], [7, 9, 18, 18], [8, 1, 30, 7], [8, 2, 25, 9], [8, 3, 5, 11], [8, 5, 20, 15], [8, 6, 30, 17], [9, 4, 12, 13], [9, 5, 24, 16], [9, 7, 18, 18]]
        node_types = [4, 5, 6, 6, 5, 2, 6, 1, 2, 3]
        node_types_ = [2, 6, 5, 1, 2, 4, 6, 3, 5, 6]
        
        # here we started indexing edges by 1, we need to change that to 0
        for edge in edge_list:
            edge[-1] = edge[-1] - 1
        
        for edge in edge_list_:
            edge[-1] = edge[-1] - 1

        code = dfs.compute_minimal_dfs_code(edge_list, node_types, timeout)
        code_ = dfs.compute_minimal_dfs_code(edge_list_, node_types_, timeout)
        self.assertTrue(len(code) == len(edge_list)//2)
        self.assertTrue(len(code_) == len(edge_list_)//2)
        self.assertTrue(isValidDFSCode(code))
        self.assertTrue(isValidDFSCode(code_))
        
        edge_dict = {r[6]: (r[5], r[7], r[3]) for r in code}  
        edge_dict_ = {r[6]: (r[5], r[7], r[3]) for r in code_}  
        for row, row_ in zip(code, code_):
		    # check dfs code
            for el, el_ in zip(row[:5], row_[:5]):
                self.assertEqual(el, el_)
		# check all vertex and edge ids
        for e, e_ in zip(edge_list, edge_list_):
            eid = e[-1]
            eid_ = e_[-1]
            de = edge_dict[eid]
            self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
            de = edge_dict_[eid_]
            e = e_
            self.assertTrue(de[2] == e[2] and ((de[0] == e[0] and de[1] == e[1]) or (de[1] == e[0] and de[0] == e[1])))
"""            
    
if __name__ == '__main__':
    unittest.main()
