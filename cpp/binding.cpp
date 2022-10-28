#include <iostream>

#include "graph_fast.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace std;
/*
 * edge_list = [(from, to, label, id)]
 */
py::list compute_minimal_dfs_code(py::list edge_list, py::list node_types, uint64_t timeout){
    Graph g;
	uint32_t id = 0;
    for (auto node_type : node_types){
        g.push_back(Vertex(node_type.cast<uint32_t>())); //liz
        //g.push_back(Vertex(node_type.cast<uint32_t>(), id)); // chris
		id++;
    }
    
    
    for (auto edge : edge_list) {
        //uint32_t from = (uint32_t) edge(0);
        vector<uint32_t> e(4);
        int idx = 0;
        for (auto el : edge) {
            e[idx] = el.cast<uint32_t>();
            //cout << e[idx] << ";";
            idx++;
        }
        g[e[0]].edges.push_back(Edge(e[0], e[1], e[2], e[3]));
    }
    g.edge_size = (int) (edge_list.size()/2);
    DFSCode code;
    get_min_code(g, code, timeout);
    std::vector<std::vector<uint32_t>> res_vec;
    
    for (auto row : code) {
        std::vector<uint32_t> vec{row.from, row.to, 
		                          row.from_label, row.e_label, row.to_label, 
		                          row.from_vertex_id, row.e_id, row.to_vertex_id};
        res_vec.push_back(vec);
    }
    
    py::list result = py::cast(res_vec);
    return result;
    
}


PYBIND11_MODULE(_dfs_codes, m) {
    m.def("compute_minimal_dfs_code", &compute_minimal_dfs_code,"Computes the minimal dfs code of a graph represented by edge_list = [(from, to, label, id)] and node_types = [type1, type2...]." 
																"Note that the entries of the list are of the form"
																 "[from_dfs_pos, to_dfs_pos,"
																 " from_label, edge_label, to_label,"
																 " from_vertex_id, edge_id, to_Vertex_id],"
																 "thus, the last three entries are not part of the DFS-code itself."
                        										 "the last three entries contain ids which allow further processing,"
																 "such as adding additional vertex and edge features to the code.");
}
