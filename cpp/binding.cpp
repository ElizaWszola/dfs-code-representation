#include <iostream>
#include "graph_fast.hpp"
#include "correlation.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace std;
/*
 * edge_list = [(from, to, label, id)]
 */
py::list compute_minimal_dfs_code(py::list edge_list, py::list node_types){
    Graph g;
    for (auto node_type : node_types){
        g.push_back(Vertex(node_type.cast<uint32_t>()));
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
    get_min_code(g, code);
    std::vector<std::vector<uint32_t>> res_vec;
    
    for (auto row : code) {
        std::vector<uint32_t> vec{row.from, row.to, row.from_label, row.e_label, row.to_label};
        res_vec.push_back(vec);
    }
    
    py::list result = py::cast(res_vec);
    return result;
    
}

py::list compute_minimal_dfs_code_vector(py::list edge_list, py::list node_types){
    Graph g;
    for (auto node_type : node_types){
        g.push_back(Vertex(node_type.cast<uint32_t>()));
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
    
    //~ std::cout << "***\n";
    //~ print_graph(g);
    //~ std::cout << "***\n";
    
    std::vector<DFSCode> codes;
    get_min_code(g, codes);
    std::vector<std::vector<uint32_t>> res_row;
    std::vector<std::vector<std::vector<uint32_t>>> res_vec;
    
    for (auto code : codes) {
      res_row.clear();
      for (auto row : code) {
          std::vector<uint32_t> vec{row.from, row.to, row.from_label, row.e_label, row.to_label};
          res_row.push_back(vec);
      }
      res_vec.push_back(res_row);
    }
    py::list result = py::cast(res_vec);
    return result;
    
}

void compute_graph_correlation(py::list node_types, py::list edge_list,
        py::array_t<double>& Y, py::array_t<bool>& B, py::array_t<bool>& fB){
    c_real* y = (c_real*)(Y.request()).ptr;
    vector<Graph> x;
    for (auto nt : node_types) {
        Graph g;
        for (auto node_type : nt){
            g.push_back(Vertex(node_type.cast<uint32_t>()));
        }
        x.push_back(g);
    }
    
    uint32_t ctr = 0;
    for (auto edl : edge_list) {
      Graph& g = x[ctr];
          uint32_t esize = 0;
          for (auto edge : edl) {
              //uint32_t from = (uint32_t) edge(0);
              vector<uint32_t> e(4);
              int idx = 0;
              for (auto el : edge) {
                  e[idx] = el.cast<uint32_t>();
                  //cout << e[idx] << ";";
                  idx++;
              }
              g[e[0]].edges.push_back(Edge(e[0], e[1], e[2], e[3]));
              ++esize;
        }
        g.edge_size = (int)(esize/2);
        ++ctr;
    }
    
    
    int64_t max_queue = 0;
    int64_t total_processed = 0;
    double delta = 0.0;
    bool* fb = (bool*)(fB.request()).ptr;
    bool* b = (bool*)(B.request()).ptr;
    
    ThreadPool tp(1);
    bool corr = tp.cmc_full(x, y, b, fb, delta, max_queue, total_processed);
}

PYBIND11_MODULE(_dfs_codes, m) {
    m.def("compute_minimal_dfs_code", &compute_minimal_dfs_code, "compute the minimal dfs code of a graph represented by edge_list = [(from, to, label, id)] and node_types = [type1, type2...]");
    m.def("compute_minimal_dfs_code_vector", &compute_minimal_dfs_code_vector, "compute the minimal dfs code of a graph represented by edge_list = [(from, to, label, id)] and node_types = [type1, type2...]");
    m.def("compute_graph_correlation", &compute_graph_correlation, "tuptuptuptuptup");
}
