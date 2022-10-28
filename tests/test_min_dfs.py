import unittest
import numpy as np
import _dfs_codes as dfs
import networkx as nx

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

class MinDFSTest(unittest.TestCase):
    
    def test_robustness_to_graph_isomorphy_nolabels(self):
        A = np.random.binomial(1, 0.3, (30, 30))
        A = np.minimum(A + A.T, 1)
        for i in range(A.shape[0]):
            A[i, i] = 0
        #node_types = np.random.randint(0, 3, A.shape[0]).tolist()
        node_types = np.ones(A.shape[0], dtype=np.int32).tolist()
        edge_list = adjacency2list(A, node_types)
        
        
        # permute the vertices and edges accordingly
        perm = np.random.permutation(A.shape[0])
        A_ = A.copy()[perm]
        A_ = A_[:, perm]
        
        g1 = nx.from_numpy_matrix(A)
        g2 = nx.from_numpy_matrix(A_)
        self.assertTrue(nx.is_isomorphic(g1, g2))
        
        node_types_ = np.asarray(node_types).copy()[perm].tolist()
        edge_list_ = adjacency2list(A_, node_types_)

        code = dfs.compute_minimal_dfs_code(edge_list, node_types)
        code_ = dfs.compute_minimal_dfs_code(edge_list_, node_types_)
        
        
        for row, row_ in zip(code, code_):
            for el, el_ in zip(row, row_):
                self.assertEqual(el, el_)
    
    
    def test_robustness_to_graph_isomorphy(self):
        
        for tries in range(100):
            A = np.random.binomial(1, 0.3, (20, 20))
            A = np.minimum(A + A.T, 1)
            for i in range(A.shape[0]):
                A[i, i] = 0
            node_types = np.random.randint(1, 7, A.shape[0]).tolist()
            edge_list = adjacency2list(A, node_types)
            
            
            # permute the vertices and edges accordingly
            perm = np.random.permutation(A.shape[0])
            A_ = A.copy()[perm]
            A_ = A_[:, perm]
            
            g1 = nx.from_numpy_matrix(A)
            g2 = nx.from_numpy_matrix(A_)
            self.assertTrue(nx.is_isomorphic(g1, g2))
            
            node_types_ = np.asarray(node_types).copy()[perm].tolist()
            edge_list_ = adjacency2list(A_, node_types_)
            


            codes = dfs.compute_minimal_dfs_code_vector(edge_list, node_types)
            codes_ = dfs.compute_minimal_dfs_code_vector(edge_list_, node_types_)
            

            code = codes[0]
            code_ = codes_[0]

            
            for row, row_ in zip(code, code_):
                for el, el_ in zip(row, row_):
                    self.assertEqual(el, el_)
                
    def test_robustness_to_graph_isomorphy_fixed(self):
        edge_list = [[0, 3, 24, 1], [0, 4, 20, 2], [0, 5, 8, 3], [0, 6, 24, 4], [0, 8, 8, 5], [1, 4, 25, 6], [1, 7, 5, 7], [2, 4, 30, 8], [2, 5, 12, 9], [2, 8, 12, 10], [2, 9, 18, 11], [3, 0, 24, 1], [3, 4, 30, 12], [3, 5, 12, 13], [4, 0, 20, 2], [4, 1, 25, 6], [4, 2, 30, 8], [4, 3, 30, 12], [4, 7, 5, 14], [5, 0, 8, 3], [5, 2, 12, 9], [5, 3, 12, 13], [5, 8, 4, 15], [6, 0, 24, 4], [6, 8, 12, 16], [6, 9, 18, 17], [7, 1, 5, 7], [7, 4, 5, 14], [7, 9, 3, 18], [8, 0, 8, 5], [8, 2, 12, 10], [8, 5, 4, 15], [8, 6, 12, 16], [9, 2, 18, 11], [9, 6, 18, 17], [9, 7, 3, 18]]
        edge_list_ = [[0, 1, 12, 1], [0, 4, 4, 2], [0, 5, 8, 3], [0, 6, 12, 4], [1, 0, 12, 1], [1, 4, 12, 5], [1, 7, 18, 6], [1, 8, 30, 7], [2, 3, 5, 8], [2, 8, 25, 9], [3, 2, 5, 8], [3, 7, 3, 10], [3, 8, 5, 11], [4, 0, 4, 2], [4, 1, 12, 5], [4, 5, 8, 12], [4, 9, 12, 13], [5, 0, 8, 3], [5, 4, 8, 12], [5, 6, 24, 14], [5, 8, 20, 15], [5, 9, 24, 16], [6, 0, 12, 4], [6, 5, 24, 14], [6, 8, 30, 17], [7, 1, 18, 6], [7, 3, 3, 10], [7, 9, 18, 18], [8, 1, 30, 7], [8, 2, 25, 9], [8, 3, 5, 11], [8, 5, 20, 15], [8, 6, 30, 17], [9, 4, 12, 13], [9, 5, 24, 16], [9, 7, 18, 18]]
        node_types = [4, 5, 6, 6, 5, 2, 6, 1, 2, 3]
        node_types_ = [2, 6, 5, 1, 2, 4, 6, 3, 5, 6]
        
        # here we started indexing edges by 1, we need to change that to 0
        for edge in edge_list:
            edge[-1] = edge[-1] - 1
        
        for edge in edge_list_:
            edge[-1] = edge[-1] - 1

        codes = dfs.compute_minimal_dfs_code_vector(edge_list, node_types)
        codes_ = dfs.compute_minimal_dfs_code_vector(edge_list_, node_types_)
        
        code = codes[0]
        code_ = codes_[0]

        
        for row, row_ in zip(code, code_):
            for el, el_ in zip(row, row_):
                self.assertEqual(el, el_)
                
    """def test_list_of_codes(self):
        # ~ return
        print("test list of codes")
        A = np.random.binomial(1, 0.3, (5, 5))
        A = np.minimum(A + A.T, 1)
        for i in range(A.shape[0]):
            A[i, i] = 0
        node_types = np.random.randint(1, 7, A.shape[0]).tolist()
        edge_list = adjacency2list(A, node_types)
        print("I run vectors and stuff")
        codes = dfs.compute_minimal_dfs_code_vector(edge_list, node_types)
        print("Done")
        for code in codes:
            print("==CODE==")
            for row in code:
              print(row)
        print("Found ", len(codes), " code candidates total.")
    """           
                
    
if __name__ == '__main__':
    unittest.main()
