import numpy as np
import _dfs_codes as dfs
import networkx as nx
import random

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


def generate_graphs_and_run():
        random.seed(1337)
        x = []
        y = []
        node_types = []
        edge_list = []
        for index in range(2):
            A = np.random.binomial(1, 0.3, (25, 25))
            A = np.minimum(A + A.T, 1)
            for i in range(A.shape[0]):
                A[i, i] = 0
            g = nx.from_numpy_matrix(A)
            if(nx.is_connected(g)):
                x.append(g)
                y.append(random.uniform(0, 1) * 100)
                node_types.append(np.random.randint(1, 7, A.shape[0]).tolist())
                edge_list.append(adjacency2list(A, node_types[-1]))
        freq = np.zeros(len(y), dtype=bool)
        fB = np.zeros(len(y), dtype=bool)
        dfs.compute_graph_correlation(node_types, edge_list, np.asarray(y), freq, fB)
                

                    
if __name__ == '__main__':
    generate_graphs_and_run()
