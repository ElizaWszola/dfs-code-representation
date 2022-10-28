# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 13:31:09 2021

@author: chris
"""

import numpy as np
import _dfs_codes as dfs
from .graph import Graph


def edgeindex_2_lists(edge_index, vertex_labels, edge_labels):
    """
    Convert pytorch geometric graph into the graph representation used for dfs code generation.
    
    Note: we assume to work with undirected graphs for now. Note that it is important that
    we feed the same edge id for both directions of the undirected edge.

    Parameters
    ----------
    data : torch_geometric.data.data.Data data
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels

    Returns
    -------
    edge_list : edge_list : [[from_vertex_id, to_vertex_id, edge_label, edge_id], ...]
    """
    edges_coo = edge_index.T
    idx_dict = {}
    idx = 0
    edge_list = []
    for eidx, (e_label, e) in enumerate(zip(edge_labels, edges_coo)):
        from_idx = e[0] 
        to_idx = e[1]
        
        key = tuple([e[1], e[0]])
        if key in idx_dict:
            e_idx = idx_dict[key]
        else:
            idx_dict[(e[0], e[1])] = idx
            e_idx = idx
            idx += 1
        edge_list += [[from_idx, to_idx, e_label, e_idx]]
    return edge_list

def min_dfs_code_from_edgeindex(edge_index, vertex_labels, edge_labels, timeout=3600):
    """
    Compute the minimal DFS code of a pytorch geometric graph.    

    Parameters
    ----------
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels
    timeout: integer timeout in seconds

    Returns
    -------
    code : minimal dfs code of the form 
    [[from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_number, to_vertex_id], ...]
    dfs_indices: [dfs_index_of_vertex_0, ...]
    """
    edges_coo = edge_index.T
    edge2id = {tuple(e.tolist()): idx for idx, e in enumerate(edges_coo)}
    edge_list = edgeindex_2_lists(edge_index, vertex_labels, edge_labels)
    code = dfs.compute_minimal_dfs_code(edge_list, vertex_labels, int(timeout*1e6))
    if code is None:
        raise TimeoutError('Timed out with a timelimit of %d seconds.'%(timeout))
    dfs_indices = {}
    for idx, row in enumerate(code):
        code[idx][-2] = edge2id[(row[-3], row[-1])]
        dfs_indices[row[-3]] = row[0]
        dfs_indices[row[-1]] = row[1]
    #print(code)
    #print(dfs_indices)
    dfs_indices = [dfs_indices[idx] for idx in range(len(vertex_labels))]
    
    return code, dfs_indices


def rnd_dfs_code_from_edgeindex(edge_index, vertex_labels, edge_labels):
    """
    Compute just any DFS code for the given graph.

    Parameters
    ----------
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels

    Returns
    -------
    code: just any dfs code of the form 
    [[from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_number, to_vertex_id], ...]
    dfs_indices: [dfs_index_of_vertex_0, ...]
    """
    edge_list = edge_index.T.tolist()
    g = Graph(edge_list, vertex_labels, edge_labels)
    return g.DFSCode()

def min_dfs_only_code_from_edgeindex(edge_index, vertex_labels, edge_labels, timeout=3600):
    """
    Compute the minimal DFS code of a pytorch geometric graph.    

    Parameters
    ----------
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels
    timeout: integer timeout in seconds

    Returns
    -------
    code : minimal dfs code of the form 
    [[from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_number, to_vertex_id], ...]
    dfs_indices: [dfs_index_of_vertex_0, ...]
    """
    edges_coo = edge_index.T
    edge2id = {tuple(e.tolist()): idx for idx, e in enumerate(edges_coo)}
    edge_list = edgeindex_2_lists(edge_index, vertex_labels, edge_labels)
    code = dfs.compute_minimal_dfs_code(edge_list, vertex_labels, int(timeout*1e6))
    if code is None:
        raise TimeoutError('Timed out with a timelimit of %d seconds.'%(timeout))
    return code


def torch_geometric_2_lists(data, vertex_labels, edge_labels):
    """
    Convert pytorch geometric graph into the graph representation used for dfs code generation.
    
    Note: we assume to work with undirected graphs for now. Note that it is important that
    we feed the same edge id for both directions of the undirected edge.

    Parameters
    ----------
    data : torch_geometric.data.data.Data data
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels

    Returns
    -------
    edge_list : edge_list : [[from_vertex_id, to_vertex_id, edge_label, edge_id], ...]
    """
    edges_coo = data.edge_index.detach().cpu().numpy().T
    idx_dict = {}
    idx = 0
    edge_list = []
    for eidx, (e_label, e) in enumerate(zip(edge_labels, edges_coo)):
        from_idx = e[0] 
        to_idx = e[1]
        
        key = tuple([e[1], e[0]])
        if key in idx_dict:
            e_idx = idx_dict[key]
        else:
            idx_dict[(e[0], e[1])] = idx
            e_idx = idx
            idx += 1
        edge_list += [[from_idx, to_idx, e_label, e_idx]]
    return edge_list


def min_dfs_code_from_torch_geometric_deprecated(data, vertex_labels, edge_labels):
    """
    Compute the minimal DFS code of a pytorch geometric graph.    

    Parameters
    ----------
    data : torch_geometric.data.data.Data data
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels

    Returns
    -------
    code : minimal dfs code of the form 
    [[from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_id, to_vertex_id], ...]
    newid2oldids: a dictionary to convert from the edge_ids required for the dfs codes to the torch geometric edge ids

    """
    edges_coo = data.edge_index.detach().cpu().numpy().T
    edge2id = {tuple(e.tolist()): idx for idx, e in enumerate(edges_coo)}
    edge_list = torch_geometric_2_lists(data, vertex_labels, edge_labels)
    code = dfs.compute_minimal_dfs_code(edge_list, vertex_labels)
    newid2oldids = {row[-2]:edge2id[(row[-3], row[-1])] for row in code}    
    return code, newid2oldids 


def min_dfs_code_from_torch_geometric_new(data, vertex_labels, edge_labels):
    return min_dfs_code_from_torch_geometric(data, vertex_labels, edge_labels)


def min_dfs_code_from_torch_geometric(data, vertex_labels, edge_labels, timeout=3600):
    """
    Compute the minimal DFS code of a pytorch geometric graph.    

    Parameters
    ----------
    data : torch_geometric.data.data.Data data
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels
    timeout: integer timeout in seconds

    Returns
    -------
    code : minimal dfs code of the form 
    [[from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_number, to_vertex_id], ...]
    dfs_indices: [dfs_index_of_vertex_0, ...]
    """
    edges_coo = data.edge_index.detach().cpu().numpy().T
    edge2id = {tuple(e.tolist()): idx for idx, e in enumerate(edges_coo)}
    edge_list = torch_geometric_2_lists(data, vertex_labels, edge_labels)
    code = dfs.compute_minimal_dfs_code(edge_list, vertex_labels, int(timeout*1e6))
    if code is None:
        raise TimeoutError('Timed out with a timelimit of %d seconds.'%(timeout))
    dfs_indices = {}
    for idx, row in enumerate(code):
        code[idx][-2] = edge2id[(row[-3], row[-1])]
        dfs_indices[row[-3]] = row[0]
        dfs_indices[row[-1]] = row[1]
    dfs_indices = [dfs_indices[idx] for idx in range(len(vertex_labels))]
    
    return code, dfs_indices


def rnd_dfs_code_from_torch_geometric(data, vertex_labels, edge_labels):
    """
    Compute just any DFS code for the given graph.

    Parameters
    ----------
    data : torch_geometric.data.data.Data data
    vertex_labels : integer list of vertex labels 
    edge_labels : integer list of edge labels

    Returns
    -------
    code: just any dfs code of the form 
    [[from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_number, to_vertex_id], ...]
    dfs_indices: [dfs_index_of_vertex_0, ...]
    """
    edge_list = data.edge_index.detach().cpu().numpy().T.tolist()
    g = Graph(edge_list, vertex_labels, edge_labels)
    return g.DFSCode()


def adjacency2list(A, node_labels):
    """
    Convert graph representation into the graph representation used for dfs code generation.

    Parameters
    ----------
    A : adjacency matrix of an undirected graph
    node_labels : list of node labels

    Returns
    -------
    edge_list : [[from_vertex_id, to_vertex_id, edge_label, edge_id], ...] with 
    edge_label = node_label[from_vertex_id] * node_label[to_vertex_id]
    """
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


def min_dfs_code(A, node_labels, timeout=3600):
    """
    Compute the minimal DFS code of a undirected vertex-labelled graph.    

    Parameters
    ----------
    A : adjacency matrix of an undirected graph
    node_labels : list of node labels
    timeout: integer timeout in seconds

    Returns
    -------
    code : minimal dfs code of the form 
    [[from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_id, to_vertex_id], ...]

    """
    edge_list = adjacency2list(A, node_labels)
    code = dfs.compute_minimal_dfs_code(edge_list, node_labels, int(timeout*1e6))
    return code



