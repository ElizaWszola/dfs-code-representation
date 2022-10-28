#ifndef GRAPH_FAST_HPP
#define GRAPH_FAST_HPP

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cassert>


#include "common.hpp"

#define MAP false

class Edge {
  public:
  uint32_t from;
  uint32_t to;
  uint32_t elabel;
  uint32_t id;
  Edge(): from(0), to(0), elabel(0), id(0) {}
  Edge(uint32_t f, uint32_t t, uint32_t l, uint32_t i):
       from(f), to(t), elabel(l), id(i) {}
};

typedef std::vector<Edge*> EdgeList;

class Vertex {
  public:
  typedef std::vector<Edge>::iterator edge_it;
  std::vector<Edge> edges;
  uint32_t label;
  Vertex(uint32_t l): label(l) {}
};

class Graph: public std::vector<Vertex> {
  public:
  uint32_t edge_size;
  Graph(): edge_size(0) {}
};

class DFSRow {
  public:
  uint32_t from;
  uint32_t to;
  uint32_t from_label;
  uint32_t e_label;
  uint32_t to_label;
  // ids
  uint32_t from_vertex_id;
  uint32_t e_id;
  uint32_t to_vertex_id;
  DFSRow(uint32_t f, uint32_t t, uint32_t fl, uint32_t el, uint32_t tl,
         uint32_t fi, uint32_t ei, uint32_t ti):
         from(f), to(t), from_label(fl), e_label(el), to_label(tl),
         from_vertex_id(fi), e_id(ei), to_vertex_id(ti) {}
  bool operator<(const DFSRow& other) {
    return
        to < from && other.to > other.from ||
        to < from && other.to < other.from && to < other.to ||
        to < from && other.to < other.from && to == other.to &&
            e_label < other.e_label ||
        to > from && other.to > other.from && other.from < from ||
        to > from && other.to > other.from && from == other.from &&
            from_label < other.from_label ||
        to > from && other.to > other.from && from == other.from &&
            from_label == other.from_label && e_label < other.e_label ||
        to > from && other.to > other.from && from == other.from &&
            from_label == other.from_label &&
            e_label == other.e_label && to_label < other.to_label;
  }
  
  bool operator==(const DFSRow& other) {
    return from == other.from && to == other.to &&
           from_label == other.from_label && e_label == other.e_label &&
           to_label == other.to_label;
  }
  bool operator!=(const DFSRow& other) {
    return !(*this == other);
  }
};

class DFSCode: public std::vector<DFSRow> {
	public:
  bool operator <(const DFSCode& other) {
     uint32_t s = size();
     uint32_t o = other.size();
     uint32_t min_size = std::min(s, o);
     for (int32_t i = 0; i < min_size; ++i) {
       if ((*this)[i] != other[i]) {
         return (*this)[i] < other[i];
       }
     }
     return o > s;
  }
  bool operator==(const DFSCode& other) {
    uint32_t s = size();
    if (s != other.size()) {
      return false;
    }
    for (int32_t i = 0; i < s; ++i) {
      if ((*this)[i] != other[i]) {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const DFSCode& other) {
    return !(*this == other);
  }
  bool is_sub(const DFSCode& other) {
    uint32_t s = size();
    uint32_t o = other.size();
    if (s > o) {
      return false;
    }
    for (int32_t i = 0; i < s; ++i) {
      if ((*this)[i] != other[i]) {
        return false;
      }
    }
    return true;
  }
};

class PDFS {
  public:
  uint32_t id;
  Edge* edge;
  PDFS* prev;
  PDFS(): id(0), edge(0), prev(0) {}	
  PDFS(uint32_t i, Edge* e, PDFS* p): id(i), edge(e), prev(p) {}	
};

class Proj: public std::vector<PDFS> {
  public:
};

class History: public std::vector<Edge*> {
  public:
  std::vector<bool> e_ids;
  std::vector<bool> v_ids;
  History(uint32_t vsize, uint32_t esize) {
    resize(esize);
    e_ids.resize(esize);
    v_ids.resize(vsize);
  }
};

typedef std::map<uint32_t, Proj> Pmap1;
typedef std::map<uint32_t, std::map<uint32_t, Proj>> Pmap2;
typedef std::map<uint32_t, std::map<
        uint32_t, std::map<uint32_t, Proj>>> Pmap3;

uint32_t get_forward_rmpath(Graph& g, Edge* e, uint32_t min_label,
                        History& history, EdgeList& edges);
uint32_t get_forward_pure(Graph& g, Edge* e, uint32_t min_label,
                      History& history, EdgeList& edges);
uint32_t get_forward_root(Graph& g, Vertex& v, EdgeList& edges);
uint32_t build_history(Graph& g, PDFS* proj, History& history);
uint32_t build_rm_path(DFSCode& code, std::vector<uint32_t>& rmpath);
void get_min_code(Graph& g, DFSCode& min_dfs);
void print_dfs_code(DFSCode& code);
void create_sample(Graph& graph);
void print_graph(Graph& graph);

void get_min_code(Graph& g, std::vector<DFSCode>& res_dfs, uint64_t timeout, bool& timed_out);
void get_min_code(Graph& g, std::vector<DFSCode>& res_dfs, uint64_t timeout);
void get_min_code(Graph& g, std::vector<DFSCode>& res_dfs);
void get_min_code(Graph& g, DFSCode& res_dfs, uint64_t timeout, bool& timed_out);
void get_min_code(Graph& g, DFSCode& res_dfs, uint64_t timeout);
void get_min_code(Graph& g, DFSCode& res_dfs);

#endif
