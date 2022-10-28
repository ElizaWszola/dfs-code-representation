#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "graph_fast.hpp"

// Get edges that start from e's "from" vertex that fulfill the
// following conditions.
// The edge is different than e (has a different "to" vertex).
// The "to" label of the edge is not smaller than min_label.
// The edge has no smaller "from edge to" label than that of e.
// The "to" is not already in history (I guess this prevents cycles).
uint32_t get_forward_rmpath(Graph& g, Edge* e, uint32_t min_label,
                        History& history, EdgeList& edges) {
  uint32_t edge_ctr = 0;
  uint32_t to_label = g[e->to].label;
  for (Vertex::edge_it i = g[e->from].edges.begin();
       i != g[e->from].edges.end(); ++i) {
    uint32_t other_to_label = g[i->to].label;
    if (e->to != i->to && min_label <= other_to_label &&
        !history.v_ids[i->to] && (e->elabel < i->elabel ||
          (e->elabel == i->elabel && to_label <= other_to_label))) {
      edges[edge_ctr] = &(*i);
      ++edge_ctr;
    }
  }
  return edge_ctr;
}

// Get edges that start from e's "to" vertex that fulfill the
// following conditions.
// The "to" label of the edge is not smaller than min_label.
// The "to" is not already in history (I guess this prevents cycles).
uint32_t get_forward_pure(Graph& g, Edge* e, uint32_t min_label,
                      History& history, EdgeList& edges) {
  uint32_t edge_ctr = 0;
  for (Vertex::edge_it i = g[e->to].edges.begin();
       i != g[e->to].edges.end(); ++i) {
    if (min_label <= g[i->to].label && !history.v_ids[i->to]) {
      edges[edge_ctr] = &(*i);
      ++edge_ctr;
    }
  }
  return edge_ctr;
}

// Get all edges that go to from the vertex v to vertices with
// less/equal labels.
uint32_t get_forward_root(Graph& g, Vertex& v, EdgeList& edges) {
  uint32_t edge_ctr = 0;
  for (Vertex::edge_it i = v.edges.begin(); i != v.edges.end(); ++i) {
    if (v.label <= g[i->to].label) {
      edges[edge_ctr] = &(*i);
      ++edge_ctr;
    }
  }
  return edge_ctr;
}

// Get one edge. It must not be in the history, it must be between
// e2 and e1, the label of e1 must be smaller than that of this edge,
// or or the label of "to" node of e1 must be no greater than the label
// of that edge's "from" node.
// Pick the first edge from between e2 and e1 that matches
// (should be unique or non-existent).
Edge* get_backward(Graph& g, Edge* e1, Edge* e2, History& history) {
  if (e1 == e2) {
    return nullptr;
  }
  for (Vertex::edge_it i = g[e2->to].edges.begin();
       i != g[e2->to].edges.end(); ++i) {
    if (!history.e_ids[i->id] && i->to == e1->from  &&
        (e1->elabel < i->elabel ||
          (e1->elabel == i->elabel &&
            g[e1->to].label <= g[e2->to].label))) {
      return &(*i);
    }
  }
  return nullptr;
}


// TODO (optimizations): rmpath looks related to history, why not use
// the same code structure in both functions below?

// Look at the projected edge e and save all its predecessors to
// history, starting from the earliest predecessor, and ending with e.
//
// TODO (optimizations): access locality screams for making class
// History an AoS. I also wonder if recursion wouldn't be better than
// iterative push_backs followed by an std::reverse.
uint32_t build_history(Graph& g, PDFS* proj, History& history) {
  uint32_t gedge_size = g.edge_size;
  uint32_t gsize = g.size();
  uint32_t hsize = 0;
  for (int32_t i = 0; i < gedge_size; ++i) {
    history.e_ids[i] = false;
  }
  for (int32_t i = 0; i < gsize; ++i) {
    history.v_ids[i] = false;
  }
  if (proj) {
    history[hsize] = proj->edge;
    ++hsize;
    history.e_ids[proj->edge->id] = true;
    history.v_ids[proj->edge->from] = true;
    history.v_ids[proj->edge->to] = true;
    for (PDFS* p = proj->prev; p; p = p->prev) {
      history[hsize] = p->edge;
      ++hsize;
      history.e_ids[p->edge->id] = true;
      history.v_ids[p->edge->from] = true;
      history.v_ids[p->edge->to] = true;
    }
    std::reverse(history.begin(), history.begin() + hsize); // TODO do we need this?
  }
  return hsize;
}

// Starting from the lowest row of the DFS code, build a path of node
// ids up to the point you can no longer go up,
// where the nodes actually form a path in the graph.
// I guess the name means "recursive mining path", but I dunno really.
uint32_t build_rm_path(DFSCode& code, std::vector<uint32_t>& rmpath) {
  uint32_t rmp_size = 0;
  // I'm using max, and not a -1, because unsigned
  uint32_t old_from = std::numeric_limits<uint32_t>::max();
  for (int32_t i = code.size() - 1; i >= 0; --i) {
    if (code[i].from < code[i].to &&
        (!rmp_size || old_from == code[i].to)) {
      rmpath[rmp_size] = i;
      ++rmp_size;
      old_from = code[i].from;
    }
  }
  return rmp_size;
}

void print_dfs_code(DFSCode& code) {
  for (int32_t i = 0; i < code.size(); ++i) {
    DFSRow row = code[i];
    std::cout << "(" << row.from << ", " << row.to << ", "
              << row.from_label << ", " << row.e_label <<  ", "
              << row.to_label << ")\n";
  }
}

// A hardcoded example graph from "Efficiently mining delta-tolerance
// closed frequent subgraphs", Fig. 2 (a).
// If the graph is not a directed one, all its edges are stored twice,
// but each pair of edges shares the id. Edge_size is the number of
// different ids.
void create_sample(Graph& graph) {
  //~ graph.push_back(Vertex(10));
  //~ graph.push_back(Vertex(10));
  //~ graph.push_back(Vertex(20));
  //~ graph.push_back(Vertex(10));
  //~ Edge er1(0, 1, 100, 0);
  //~ Edge el1(1, 0, 100, 0);
  //~ Edge er2(0, 2, 100, 1);
  //~ Edge el2(2, 0, 100, 1);
  //~ Edge er3(1, 2, 100, 2);
  //~ Edge el3(2, 1, 100, 2);
  //~ Edge er4(2, 3, 100, 3);
  //~ Edge el4(3, 2, 100, 3);
  //~ graph[0].edges.push_back(er1);
  //~ graph[0].edges.push_back(er2);
  //~ graph[1].edges.push_back(el1);
  //~ graph[1].edges.push_back(er3);
  //~ graph[2].edges.push_back(el2);
  //~ graph[2].edges.push_back(el3);
  //~ graph[2].edges.push_back(er4);
  //~ graph[3].edges.push_back(el4);
  //~ graph.edge_size = 4;
  graph.push_back(Vertex(4));
  graph.push_back(Vertex(5));
  graph.push_back(Vertex(6));
  graph.push_back(Vertex(6));
  graph.push_back(Vertex(5));
  graph.push_back(Vertex(2));
  graph.push_back(Vertex(6));
  graph.push_back(Vertex(1));
  graph.push_back(Vertex(2));
  graph.push_back(Vertex(3));
  std::vector<Edge> edges;
  {
    Edge e(0, 3, 24, 1);
    edges.push_back(e);
  } {
    Edge e(0, 4, 20, 2);
    edges.push_back(e);
  } {
  Edge e(0, 5, 8, 3);
    edges.push_back(e);
  } {
  Edge e(0, 6, 24, 4);
    edges.push_back(e);
  } {
  Edge e(0, 8, 8, 5);
    edges.push_back(e);
  } {
  Edge e(1, 4, 25, 6);
    edges.push_back(e);
  } {
  Edge e(1, 7, 5, 7);
    edges.push_back(e);
  } {
  Edge e(2, 4, 30, 8);
    edges.push_back(e);
  } {
  Edge e(2, 5, 12, 9);
    edges.push_back(e);
  } {
  Edge e(2, 8, 12, 10);
    edges.push_back(e);
  } {
  Edge e(2, 9, 18, 11);
    edges.push_back(e);
  } {
  Edge e(3, 0, 24, 1);
    edges.push_back(e);
  } {
  Edge e(3, 4, 30, 12);
    edges.push_back(e);
  } {
  Edge e(3, 5, 12, 13);
    edges.push_back(e);
  } {
  Edge e(4, 0, 20, 2);
    edges.push_back(e);
  } {
  Edge e(4, 1, 25, 6);
    edges.push_back(e);
  } {
  Edge e(4, 2, 30, 8);
    edges.push_back(e);
  } {
  Edge e(4, 3, 30, 12);
    edges.push_back(e);
  } {
  Edge e(4, 7, 5, 14);
    edges.push_back(e);
  } {
  Edge e(5, 0, 8, 3);
    edges.push_back(e);
  } {
  Edge e(5, 2, 12, 9);
    edges.push_back(e);
  } {
  Edge e(5, 3, 12, 13);
    edges.push_back(e);
  } {
  Edge e(5, 8, 4, 15);
    edges.push_back(e);
  } {
  Edge e(6, 0, 24, 4);
    edges.push_back(e);
  } {
  Edge e(6, 8, 12, 16);
    edges.push_back(e);
  } {
  Edge e(6, 9, 18, 17);
    edges.push_back(e);
  } {
  Edge e(7, 1, 5, 7);
    edges.push_back(e);
  } {
  Edge e(7, 4, 5, 14);
    edges.push_back(e);
  } {
  Edge e(7, 9, 3, 18);
    edges.push_back(e);
  } {
  Edge e(8, 0, 8, 5);
    edges.push_back(e);
  } {
  Edge e(8, 2, 12, 10);
    edges.push_back(e);
  } {
  Edge e(8, 5, 4, 15);
    edges.push_back(e);
  } {
  Edge e(8, 6, 12, 16);
    edges.push_back(e);
  } {
  Edge e(9, 2, 18, 11);
    edges.push_back(e);
  } {
  Edge e(9, 6, 18, 17);
    edges.push_back(e);
  } {
  Edge e(9, 7, 3, 18);
    edges.push_back(e);
  }
  for (int32_t i = 0; i < edges.size(); ++i) {
    Edge e = edges[i];
    e.id = e.id - 1;
    graph[e.from].edges.push_back(e);
  }
  graph.edge_size = 18;
}

// Just a helper function to see if what I made makes sense.
void print_graph(Graph& graph) {
  for (int32_t i = 0; i < graph.size(); ++i) {
    std::cout << i << " " << graph[i].label << "\n";
    for (int32_t j = 0; j < graph[i].edges.size(); ++j) {
      Edge& e = graph[i].edges[j];
      std::cout << "\t" << e.id << " " << e.from << " "
                << e.to << " " << e.elabel << "\n";
    }
  }
}

inline bool is_label_eq(uint32_t from1, uint32_t edge1, uint32_t to1,
                        uint32_t from2, uint32_t edge2, uint32_t to2) {
  return from1 == from2 && edge1 == edge2 && to1 == to2;
}

inline bool is_label_lt(uint32_t from1, uint32_t edge1, uint32_t to1,
                        uint32_t from2, uint32_t edge2, uint32_t to2) {
  return from1 < from2 ||
         from1 == from2 && (edge1 < edge2 ||
                            edge1 == edge2 && to1 < to2);
}

inline bool is_label_eq(uint32_t edge1, uint32_t to1,
                        uint32_t edge2, uint32_t to2) {
  return edge1 == edge2 && to1 == to2;
}

inline bool is_label_lt(uint32_t edge1, uint32_t to1,
                        uint32_t edge2, uint32_t to2) {
  return edge1 < edge2 || edge1 == edge2 && to1 < to2;
}

inline bool is_label_eq(uint32_t edge1, uint32_t edge2) {
  return edge1 == edge2;
}

inline bool is_label_lt(uint32_t edge1, uint32_t edge2) {
  return edge1 < edge2;
}


// Min_dfs actually stands for a min DFS candidate.
// To my understanding, this function simply gradually reconstructs g
// from scratch in the order defined by min dfs. The added edges are
// saved in the creation order in vector projected.
// The order of building is like so:
// 1. Try to reconstruct backward connections starting from the least
//    recently added node up to the most recently added.
// 2. Try to reconstruct forward connections from the most recenly
//    added node.
// 3. Try to reconstruct forward connections starting from the most
//    recently added node down to the least recently added.
//
void get_min_code_rec(Graph& g, DFSCode& min_dfs, PDFS& cur_edge,
                      std::vector<uint32_t>& rmpath, History& history,
                      EdgeList& edges, std::vector<DFSCode>& res_dfs,
                      uint64_t timeout, clocktime &start) {
  if (timeout) {
    clocktime current;
    get_time(&current);
    if (get_time_difference(&start, &current) > timeout * 1000) {
      res_dfs.clear();
      return;
    }
  }
  uint32_t rmp_size = build_rm_path(min_dfs, rmpath);
  // the "from" label of the first edge of the dfs code
  uint32_t min_label = min_dfs[0].from_label;
  // the "to" id of the last forward edge of the dfs code
  uint32_t max_id = min_dfs[rmpath[0]].to;
  #if MAP
  Pmap1 root;
  #else
  Proj projected;
  uint32_t min_edge = std::numeric_limits<uint32_t>::max();
  uint32_t min_to = std::numeric_limits<uint32_t>::max();
  #endif
  bool added = false; // added new edge?
  uint32_t new_to = 0;
  // generate backward edges
  uint32_t hsize = build_history(g, &cur_edge, history);
  for (int32_t i = rmp_size - 1; !added && i >= 1; --i) {
    // e1 is expected to be earlier in the DFSCode than e2,
    // because rmpath is built "bottom up".
    Edge *e = get_backward(g, history[rmpath[i]],
                           history[rmpath[0]], history);
    if (e) {
      #if MAP
      root[e->elabel].push_back(PDFS(0, e, &cur_edge));
      #else
      uint32_t c_edge = e->elabel;
      if (is_label_eq(c_edge, min_edge)) {
        projected.push_back(PDFS(0, e, &cur_edge));
      } else if (is_label_lt(c_edge, min_edge)) {
        min_edge = c_edge;
        projected.clear();
        projected.push_back(PDFS(0, e, &cur_edge));
      }
      #endif
      new_to = min_dfs[rmpath[i]].from;
      added = true;
    }
  }
  if (added) {
    #if MAP
    Proj& projected = root.begin()->second;
    #endif
    for (int32_t i = 1; i < projected.size(); ++i) {
      #if MAP
      DFSRow new_row = DFSRow(max_id, new_to, g[projected[i].edge->from].label,
                              root.begin()->first, g[projected[i].edge->to].label);
      #else
      DFSRow new_row = DFSRow(max_id, new_to, g[projected[i].edge->from].label,
                              min_edge, g[projected[i].edge->to].label);
      #endif
      DFSCode min_dfs_local(min_dfs);
      min_dfs_local.push_back(new_row);
      std::vector<uint32_t> remove;
      for (int32_t i = 0; i < res_dfs.size(); ++i) {
        DFSCode& other = res_dfs[i];
        if (!min_dfs_local.is_sub(other)) {
          if (min_dfs_local < other) {
            remove.push_back(i);
          } else if (other < min_dfs_local) {
            return;
          }
        }
      }
      for (int32_t i = remove.size() - 1; i >= 0; --i) {
        res_dfs.erase(res_dfs.begin() + remove[i]);
      }
      get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                       edges, res_dfs, timeout, start);
    }
    #if MAP
    DFSRow new_row = DFSRow(max_id, new_to, g[projected[0].edge->from].label,
                               root.begin()->first, g[projected[0].edge->to].label);
    #else
    DFSRow new_row = DFSRow(max_id, new_to, g[projected[0].edge->from].label,
                               min_edge, g[projected[0].edge->to].label);
    #endif
    min_dfs.push_back(new_row);
    std::vector<uint32_t> remove;
    for (int32_t i = 0; i < res_dfs.size(); ++i) {
      DFSCode& other = res_dfs[i];
      if (!min_dfs.is_sub(other)) {
        if (min_dfs < other) {
          remove.push_back(i);
        } else if (other < min_dfs) {
          return;
        }
      }
    }
    for (int32_t i = remove.size() - 1; i >= 0; --i) {
      res_dfs.erase(res_dfs.begin() + remove[i]);
    }
    get_min_code_rec(g, min_dfs, projected[0], rmpath, history,
                    edges, res_dfs, timeout, start);
    return;
  }
  #if MAP
  Pmap2 root2;
  #endif
  uint32_t new_from = 0;
  added = false;
  //generate edges from the last added edge
  uint32_t edge_size = get_forward_pure(g, history[rmpath[0]],
                                        min_label, history, edges);
  if (edge_size) {
    added = true;
    new_from = max_id;
    for (uint32_t eid = 0; eid < edge_size; ++eid) {
      Edge* e = edges[eid];
      #if MAP
      root2[e->elabel][g[e->to].label].push_back(
          PDFS(0, e, &cur_edge));
      #else
      uint32_t c_edge = e->elabel;
      uint32_t c_to = g[e->to].label;
      if (is_label_eq(c_edge, c_to, min_edge, min_to)) {
        projected.push_back(PDFS(0, e, &cur_edge));
      } else if (is_label_lt(c_edge, c_to, min_edge, min_to)) {
        min_edge = c_edge;
        min_to = c_to;
        projected.clear();
        projected.push_back(PDFS(0, e, &cur_edge));
      }
      #endif
      new_from = max_id;
    }
  }
  // generate forward edges
  for (int32_t i = 0; !added && i < rmp_size; ++i) {
    uint32_t edge_size = get_forward_rmpath(g, history[rmpath[i]],
                                            min_label, history, edges);
    if (edge_size) {
      added = true;
      new_from = min_dfs[rmpath[i]].from;
      for (uint32_t eid = 0; eid < edge_size; ++eid) {
        Edge* e = edges[eid];
        #if MAP
        root2[e->elabel][g[e->to].label].push_back(
            PDFS(0, e, &cur_edge));
        #else
        uint32_t c_edge = e->elabel;
        uint32_t c_to = g[e->to].label;
        if (is_label_eq(c_edge, c_to, min_edge, min_to)) {
          projected.push_back(PDFS(0, e, &cur_edge));
        } else if (is_label_lt(c_edge, c_to, min_edge, min_to)) {
          min_edge = c_edge;
          min_to = c_to;
          projected.clear();
          projected.push_back(PDFS(0, e, &cur_edge));
        }
        #endif
        new_from = min_dfs[rmpath[i]].from;
      }
    }
  }
  if (added) {
    #if MAP
    Proj& projected = root2.begin()->second.begin()->second;
    #endif
    for (int32_t i = 1; i < projected.size(); ++i) {
      #if MAP
      DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[i].edge->from].label,
                              root2.begin()->first,
                              root2.begin()->second.begin()->first);
      #else
      DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[i].edge->from].label,
                              min_edge, min_to);
      #endif
      DFSCode min_dfs_local(min_dfs);
      min_dfs_local.push_back(new_row);
      std::vector<uint32_t> remove;
      for (int32_t i = 0; i < res_dfs.size(); ++i) {
        DFSCode& other = res_dfs[i];
        if (!min_dfs_local.is_sub(other)) {
          if (min_dfs_local < other) {
            remove.push_back(i);
          } else if (other < min_dfs_local) {
            return;
          }
        }
      }
      for (int32_t i = remove.size() - 1; i >= 0; --i) {
        res_dfs.erase(res_dfs.begin() + remove[i]);
      }
      get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                       edges, res_dfs, timeout, start);
    }
    #if MAP
    DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[0].edge->from].label,
                             root2.begin()->first,
                             root2.begin()->second.begin()->first);
    #else
    DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[0].edge->from].label,
                            min_edge, min_to);
    #endif
    min_dfs.push_back(new_row);
    std::vector<uint32_t> remove;
    for (int32_t i = 0; i < res_dfs.size(); ++i) {
      DFSCode& other = res_dfs[i];
      if (!min_dfs.is_sub(other)) {
        if (min_dfs < other) {
          remove.push_back(i);
        } else if (other < min_dfs) {
          return;
        }
      }
    }
    for (int32_t i = remove.size() - 1; i >= 0; --i) {
      res_dfs.erase(res_dfs.begin() + remove[i]);
    }
    get_min_code_rec(g, min_dfs, projected[0], rmpath, history,
                     edges, res_dfs, timeout, start);
    return;
  }
  //~ std::cout << "Code candidate:\n";
  //~ print_dfs_code(min_dfs);
  res_dfs.push_back(min_dfs);
}

// Assign to each "from-edge-to" label triplet a list of all 
// corresponding forward edges in the graph, add the first row from the
// first "from-edge-to" triplet to your dfs code and proceed recursively.
// This basically generates all possible forward edges and picks one
// with the smallest "from-edge-to" label.
//
// TODO (optimization): no need to save all edges, it might be good
// to implement a custom label compare function and always save just the
// smallest one.
// TODO what about multiple identical min labels? This code ignores this
// possibility...
void get_min_code(Graph& g, std::vector<DFSCode>& res_dfs, uint64_t timeout) {
  clocktime start;
  get_time(&start);
  EdgeList edges(g.edge_size);
  std::vector<uint32_t> rmpath(g.edge_size);
  History history(g.size(), g.edge_size);
  #if MAP
  Pmap3 root;
  #else
  Proj projected;
  uint32_t min_from = std::numeric_limits<uint32_t>::max();
  uint32_t min_edge = std::numeric_limits<uint32_t>::max();
  uint32_t min_to = std::numeric_limits<uint32_t>::max();
  #endif
  for (int32_t i = 0; i < g.size(); ++i) {
    uint32_t edge_ctr = get_forward_root(g, g[i], edges);
    for (uint32_t eid = 0; eid < edge_ctr; ++eid) {
      Edge* e = edges[eid];
      #if MAP
      root[g[i].label][e->elabel][g[e->to].label].push_back(
          PDFS(0, e, nullptr));
      #else
      uint32_t c_from = g[i].label;
      uint32_t c_edge = e->elabel;
      uint32_t c_to = g[e->to].label;
      if (is_label_eq(c_from, c_edge, c_to,
                      min_from, min_edge, min_to)) {
        projected.push_back(PDFS(0, e, nullptr));
      } else if (is_label_lt(c_from, c_edge, c_to,
                             min_from, min_edge, min_to)) {
        min_from = c_from;
        min_edge = c_edge;
        min_to = c_to;
        projected.clear();
        projected.push_back(PDFS(0, e, nullptr));
      }
      #endif
    }
  }
  #if MAP
  Proj& projected =
      root.begin()->second.begin()->second.begin()->second;
  for (int32_t i = 0; i < projected.size(); ++i) {
    DFSCode min_dfs_local;
    min_dfs_local.reserve(g.edge_size);
    min_dfs_local.push_back(
      DFSRow(0, 1, root.begin()->first,
             root.begin()->second.begin()->first,
             root.begin()->second.begin()->second.begin()->first));
    get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                     edges, res_dfs, timeout, start);
  }
  #else
  for (int32_t i = 0; i < projected.size(); ++i) {
    DFSCode min_dfs_local;
    min_dfs_local.reserve(g.edge_size);
    min_dfs_local.push_back(DFSRow(0, 1, min_from, min_edge, min_to));
    get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                     edges, res_dfs, timeout, start);
  }
  #endif
}

// Min_dfs actually stands for a min DFS candidate.
// To my understanding, this function simply gradually reconstructs g
// from scratch in the order defined by min dfs. The added edges are
// saved in the creation order in vector projected.
// The order of building is like so:
// 1. Try to reconstruct backward connections starting from the least
//    recently added node up to the most recently added.
// 2. Try to reconstruct forward connections from the most recenly
//    added node.
// 3. Try to reconstruct forward connections starting from the most
//    recently added node down to the least recently added.
//
void get_min_code_rec(Graph& g, DFSCode& min_dfs, PDFS& cur_edge,
                      std::vector<uint32_t>& rmpath, History& history,
                      EdgeList& edges, DFSCode& res_dfs, uint64_t timeout,
                      clocktime &start) {
  if (timeout) {
    clocktime current;
    get_time(&current);
    if (get_time_difference(&start, &current) > timeout * 1000) {
      res_dfs.clear();
      return;
    }
  }
  uint32_t rmp_size = build_rm_path(min_dfs, rmpath);
  // the "from" label of the first edge of the dfs code
  uint32_t min_label = min_dfs[0].from_label;
  // the "to" id of the last forward edge of the dfs code
  uint32_t max_id = min_dfs[rmpath[0]].to;
  #if MAP
  Pmap1 root;
  #else
  Proj projected;
  uint32_t min_edge = std::numeric_limits<uint32_t>::max();
  uint32_t min_to = std::numeric_limits<uint32_t>::max();
  #endif
  bool added = false; // added new edge?
  uint32_t new_to = 0;
  // generate backward edges
  uint32_t hsize = build_history(g, &cur_edge, history);
  for (int32_t i = rmp_size - 1; !added && i >= 1; --i) {
    // e1 is expected to be earlier in the DFSCode than e2,
    // because rmpath is built "bottom up".
    Edge *e = get_backward(g, history[rmpath[i]],
                           history[rmpath[0]], history);
    if (e) {
      #if MAP
      root[e->elabel].push_back(PDFS(0, e, &cur_edge));
      #else
      uint32_t c_edge = e->elabel;
      if (is_label_eq(c_edge, min_edge)) {
        projected.push_back(PDFS(0, e, &cur_edge));
      } else if (is_label_lt(c_edge, min_edge)) {
        min_edge = c_edge;
        projected.clear();
        projected.push_back(PDFS(0, e, &cur_edge));
      }
      #endif
      new_to = min_dfs[rmpath[i]].from;
      added = true;
    }
  }
  if (added) {
    #if MAP
    Proj& projected = root.begin()->second;
    #endif
    for (int32_t i = 1; i < projected.size(); ++i) {
      #if MAP
      DFSRow new_row = DFSRow(max_id, new_to, g[projected[i].edge->from].label,
                              root.begin()->first, g[projected[i].edge->to].label);
      #else
      DFSRow new_row = DFSRow(max_id, new_to, g[projected[i].edge->from].label,
                              min_edge, g[projected[i].edge->to].label);
      #endif
      DFSCode min_dfs_local(min_dfs);
      min_dfs_local.push_back(new_row);
      if (!res_dfs.empty() && res_dfs < min_dfs_local) {
        return;
      }
      get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                       edges, res_dfs, timeout, start);
    }
    #if MAP
    DFSRow new_row = DFSRow(max_id, new_to, g[projected[0].edge->from].label,
                               root.begin()->first, g[projected[0].edge->to].label);
    #else
    DFSRow new_row = DFSRow(max_id, new_to, g[projected[0].edge->from].label,
                               min_edge, g[projected[0].edge->to].label);
    #endif
    min_dfs.push_back(new_row);
    if (!res_dfs.empty() && res_dfs < min_dfs) {
      return;
    }
    get_min_code_rec(g, min_dfs, projected[0], rmpath, history,
                    edges, res_dfs, timeout, start);
    return;
  }
  #if MAP
  Pmap2 root2;
  #endif
  uint32_t new_from = 0;
  added = false;
  //generate edges from the last added edge
  uint32_t edge_size = get_forward_pure(g, history[rmpath[0]],
                                        min_label, history, edges);
  if (edge_size) {
    added = true;
    new_from = max_id;
    for (uint32_t eid = 0; eid < edge_size; ++eid) {
      Edge* e = edges[eid];
      #if MAP
      root2[e->elabel][g[e->to].label].push_back(
          PDFS(0, e, &cur_edge));
      #else
      uint32_t c_edge = e->elabel;
      uint32_t c_to = g[e->to].label;
      if (is_label_eq(c_edge, c_to, min_edge, min_to)) {
        projected.push_back(PDFS(0, e, &cur_edge));
      } else if (is_label_lt(c_edge, c_to, min_edge, min_to)) {
        min_edge = c_edge;
        min_to = c_to;
        projected.clear();
        projected.push_back(PDFS(0, e, &cur_edge));
      }
      #endif
      new_from = max_id;
    }
  }
  // generate forward edges
  for (int32_t i = 0; !added && i < rmp_size; ++i) {
    uint32_t edge_size = get_forward_rmpath(g, history[rmpath[i]],
                                            min_label, history, edges);
    if (edge_size) {
      added = true;
      new_from = min_dfs[rmpath[i]].from;
      for (uint32_t eid = 0; eid < edge_size; ++eid) {
        Edge* e = edges[eid];
        #if MAP
        root2[e->elabel][g[e->to].label].push_back(
            PDFS(0, e, &cur_edge));
        #else
        uint32_t c_edge = e->elabel;
        uint32_t c_to = g[e->to].label;
        if (is_label_eq(c_edge, c_to, min_edge, min_to)) {
          projected.push_back(PDFS(0, e, &cur_edge));
        } else if (is_label_lt(c_edge, c_to, min_edge, min_to)) {
          min_edge = c_edge;
          min_to = c_to;
          projected.clear();
          projected.push_back(PDFS(0, e, &cur_edge));
        }
        #endif
        new_from = min_dfs[rmpath[i]].from;
      }
    }
  }
  if (added) {
    #if MAP
    Proj& projected = root2.begin()->second.begin()->second;
    #endif
    for (int32_t i = 1; i < projected.size(); ++i) {
      #if MAP
      DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[i].edge->from].label,
                              root2.begin()->first,
                              root2.begin()->second.begin()->first);
      #else
      DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[i].edge->from].label,
                              min_edge, min_to);
      #endif
      DFSCode min_dfs_local(min_dfs);
      min_dfs_local.push_back(new_row);
      if (!res_dfs.empty() && res_dfs < min_dfs_local) {
        return;
      }
      get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                       edges, res_dfs, timeout, start);
    }
    #if MAP
    DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[0].edge->from].label,
                             root2.begin()->first,
                             root2.begin()->second.begin()->first);
    #else
    DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[0].edge->from].label,
                            min_edge, min_to);
    #endif
    min_dfs.push_back(new_row);
    if (!res_dfs.empty() && res_dfs < min_dfs) {
      return;
    }
    get_min_code_rec(g, min_dfs, projected[0], rmpath, history,
                     edges, res_dfs, timeout, start);
    return;
  }
  //~ std::cout << "Code candidate:\n";
  //~ print_dfs_code(min_dfs);
  res_dfs = min_dfs;
}

// timeout in ms: if timed-out, returns an empty result
void get_min_code(Graph& g, DFSCode& res_dfs, uint64_t timeout) {
  clocktime start;
  get_time(&start);
  res_dfs.reserve(g.edge_size);
  EdgeList edges(g.edge_size);
  std::vector<uint32_t> rmpath(g.edge_size);
  History history(g.size(), g.edge_size);
  #if MAP
  Pmap3 root;
  #else
  Proj projected;
  projected.reserve(g.edge_size);
  uint32_t min_from = std::numeric_limits<uint32_t>::max();
  uint32_t min_edge = std::numeric_limits<uint32_t>::max();
  uint32_t min_to = std::numeric_limits<uint32_t>::max();
  #endif
  for (int32_t i = 0; i < g.size(); ++i) {
    uint32_t edge_ctr = get_forward_root(g, g[i], edges);
    for (uint32_t eid = 0; eid < edge_ctr; ++eid) {
      Edge* e = edges[eid];
      #if MAP
      root[g[i].label][e->elabel][g[e->to].label].push_back(
          PDFS(0, e, nullptr));
      #else
      uint32_t c_from = g[i].label;
      uint32_t c_edge = e->elabel;
      uint32_t c_to = g[e->to].label;
      if (is_label_eq(c_from, c_edge, c_to,
                      min_from, min_edge, min_to)) {
        projected.push_back(PDFS(0, e, nullptr));
      } else if (is_label_lt(c_from, c_edge, c_to,
                             min_from, min_edge, min_to)) {
        min_from = c_from;
        min_edge = c_edge;
        min_to = c_to;
        projected.clear();
        projected.push_back(PDFS(0, e, nullptr));
      }
      #endif
    }
  }
  #if MAP
  Proj& projected =
      root.begin()->second.begin()->second.begin()->second;
  for (int32_t i = 0; i < projected.size(); ++i) {
    DFSCode min_dfs_local;
    min_dfs_local.reserve(g.edge_size);
    min_dfs_local.push_back(
      DFSRow(0, 1, root.begin()->first,
             root.begin()->second.begin()->first,
             root.begin()->second.begin()->second.begin()->first));
    get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                     edges, res_dfs, timeout, start);
  }
  #else
  for (int32_t i = 0; i < projected.size(); ++i) {
    DFSCode min_dfs_local;
    min_dfs_local.reserve(g.edge_size);
    min_dfs_local.push_back(DFSRow(0, 1, min_from, min_edge, min_to));
    get_min_code_rec(g, min_dfs_local, projected[i], rmpath, history,
                     edges, res_dfs, timeout, start);
  }
  #endif
}

void get_min_code(Graph& g, DFSCode& res_dfs) {
  get_min_code(g, res_dfs, 0);
}

/*
void compute_max_correlation_rec(Graph& g, DFSCode& code, PDFS& cur_edge,
                      std::vector<uint32_t>& rmpath, History& history,
                      EdgeList& edges) {
  uint32_t rmp_size = build_rm_path(min_dfs, rmpath);
  // the "from" label of the first edge of the dfs code
  uint32_t min_label = min_dfs[0].from_label;
  // the "to" id of the last forward edge of the dfs code
  uint32_t max_id = min_dfs[rmpath[0]].to;
  //~ #if MAP
  Pmap3 fwd_root;
  Pmap2 bwd_root;
  //~ #else
  //~ Proj projected;
  //~ uint32_t min_edge = std::numeric_limits<uint32_t>::max();
  //~ uint32_t min_to = std::numeric_limits<uint32_t>::max();
  //~ #endif
  bool added = false; // added new edge?
  uint32_t new_to = 0;
  // generate backward edges
  uint32_t hsize = build_history(g, &cur_edge, history);
  for (int32_t i = rmp_size - 1; i >= 1; --i) {
    Edge *e = get_backward(g, history[rmpath[i]],
                           history[rmpath[0]], history);
    if (e) {
      //~ #if MAP
      bwd_root[g[e->from].label][e->elabel].push_back(PDFS(0, e, &cur_edge));
      //~ #else
      //~ uint32_t c_edge = e->elabel;
      //~ if (is_label_eq(c_edge, min_edge)) {
        //~ projected.push_back(PDFS(0, e, &cur_edge));
      //~ } else if (is_label_lt(c_edge, min_edge)) {
        //~ min_edge = c_edge;
        //~ projected.clear();
        //~ projected.push_back(PDFS(0, e, &cur_edge));
      //~ }
      //~ #endif
    }
  }
  //generate edges from the last added edge
  uint32_t edge_size = get_forward_pure(g, history[rmpath[0]],
                                        min_label, history, edges);
  for (uint32_t eid = 0; eid < edge_size; ++eid) {
    Edge* e = edges[eid];
    //~ #if MAP
    fwd_root[g[e->from].label][e->elabel][g[e->to].label].push_back(
        PDFS(0, e, &cur_edge));
    //~ #else
    //~ uint32_t c_edge = e->elabel;
    //~ uint32_t c_to = g[e->to].label;
    //~ if (is_label_eq(c_edge, c_to, min_edge, min_to)) {
      //~ projected.push_back(PDFS(0, e, &cur_edge));
    //~ } else if (is_label_lt(c_edge, c_to, min_edge, min_to)) {
      //~ min_edge = c_edge;
      //~ min_to = c_to;
      //~ projected.clear();
      //~ projected.push_back(PDFS(0, e, &cur_edge));
    //~ }
    //~ #endif
  }
  // generate forward edges
  for (int32_t i = 0; i < rmp_size; ++i) {
    uint32_t edge_size = get_forward_rmpath(g, history[rmpath[i]],
                                            min_label, history, edges);
    for (uint32_t eid = 0; eid < edge_size; ++eid) {
      Edge* e = edges[eid];
      //~ #if MAP
      fwd_root[g[e->from].label][e->elabel][g[e->to].label].push_back(
          PDFS(0, e, &cur_edge));
      //~ #else
      //~ uint32_t c_edge = e->elabel;
      //~ uint32_t c_to = g[e->to].label;
      //~ if (is_label_eq(c_edge, c_to, min_edge, min_to)) {
        //~ projected.push_back(PDFS(0, e, &cur_edge));
      //~ } else if (is_label_lt(c_edge, c_to, min_edge, min_to)) {
        //~ min_edge = c_edge;
        //~ min_to = c_to;
        //~ projected.clear();
        //~ projected.push_back(PDFS(0, e, &cur_edge));
      //~ }
      //~ #endif
    }
  }
  
  // generate from maps
  for (Pmap2::iterator to = bwd_root.begin(); to != bwd_root.end(); ++to) {
    for (Pmap1::iterator elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
      code.push_back(max_id, to->first, -1, elabel->first, -1);
      generate_rec(g, code, elabel->second, rmpath, history, edges);
      code.pop();
    }
  }
  for (Pmap3::iterator from = fwd_root.begin(); from != fwd_root.end(); ++from) {
    for (Pmap2::iterator elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
      for (Pmap1::iterator to = elabel->second.begin(); to != elabel->second.end(); ++to) {
        code.push_back(from->first, max_id + 1, -1, elabel->first, -1);
        generate_rec(g, code, to->second, rmpath, history, edges);
        code.pop();
      }
    }
  }
}

void generate(Graph& g, Proj& projected, DFSCode& code) {
  DFSCode mincode;
  // TODO check if code is minimal
  
}
*/
