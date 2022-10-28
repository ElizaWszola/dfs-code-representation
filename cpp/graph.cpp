#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include "graph.hpp"

// Get edges that start from e's "from" vertex that fulfill the
// following conditions.
// The edge is different than e (has a different "to" vertex).
// The "to" label of the edge is not smaller than min_label.
// The edge has no smaller "from edge to" label than that of e.
// The "to" is not already in history (I guess this prevents cycles).
bool get_forward_rmpath(Graph& g, Edge* e, uint32_t min_label,
                        History& history, EdgeList& edges) {
  edges.clear();
  uint32_t to_label = g[e->to].label;
  for (Vertex::edge_it i = g[e->from].edges.begin();
       i != g[e->from].edges.end(); ++i) {
    uint32_t other_to_label = g[i->to].label;
    if (e->to != i->to && min_label <= other_to_label &&
        !history.v_ids[i->to] && (e->elabel < i->elabel ||
          (e->elabel == i->elabel && to_label <= other_to_label))) {
      edges.push_back(&(*i));
    }
  }
  return !edges.empty();
}

// Get edges that start from e's "to" vertex that fulfill the
// following conditions.
// The "to" label of the edge is not smaller than min_label.
// The "to" is not already in history (I guess this prevents cycles).
bool get_forward_pure(Graph& g, Edge* e, uint32_t min_label,
                      History& history, EdgeList& edges) {
  edges.clear();
  for (Vertex::edge_it i = g[e->to].edges.begin();
       i != g[e->to].edges.end(); ++i) {
    if (min_label <= g[i->to].label && !history.v_ids[i->to]) {
      edges.push_back(&(*i));
    }
  }
  return !edges.empty();
}

// Get all edges that go to from the vertex v to vertices with
// less/equal labels.
bool get_forward_root(Graph& g, Vertex& v, EdgeList& edges) {
  edges.clear();
  for (Vertex::edge_it i = v.edges.begin(); i != v.edges.end(); ++i) {
    if (v.label <= g[i->to].label) {
      edges.push_back(&(*i));
    }
  }
  return !edges.empty();
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
void build_history(Graph& g, PDFS* proj, History& history) {
  history.clear();
  history.e_ids.clear();
  history.v_ids.clear();
  history.e_ids.resize(g.edge_size, 0);
  history.v_ids.resize(g.size(), 0);
  if (proj) {
    history.push_back(proj->edge); // this is apparently very expensive
    //~ assert (proj->edge->id < g.edge_size && proj->edge->from < g.size() &&
            //~ proj->edge->to < g.size());
    history.e_ids[proj->edge->id] = 1;
    history.v_ids[proj->edge->from] = 1;
    history.v_ids[proj->edge->to] = 1;
    for (PDFS* p = proj->prev; p; p = p->prev) {
      history.push_back(p->edge);
      history.e_ids[p->edge->id] = 1;
      history.v_ids[p->edge->from] = 1;
      history.v_ids[p->edge->to] = 1;
    }
    std::reverse(history.begin(), history.end());
  }
}

// Starting from the lowest row of the DFS code, build a path of node
// ids up to the point you can no longer go up,
// where the nodes actually form a path in the graph.
// I guess the name means "recursive mining path", but I dunno really.
void build_rm_path(DFSCode& code, std::vector<uint32_t>& rmpath) {
  rmpath.clear();
  // I'm using max, and not a -1, because unsigned
  uint32_t old_from = std::numeric_limits<uint32_t>::max();
  for (int32_t i = code.size() - 1; i >= 0; --i) {
    if (code[i].from < code[i].to &&
        (rmpath.empty() || old_from == code[i].to)) {
      rmpath.push_back(i);
      old_from = code[i].from;
    }
  }
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
// TODO check if the min label edge can be picked arbitrary without
// violating the minimality requirement later on. If it's not the case,
// implement the function the way you operate over multiple DFS copies
// and then pick the best one. Would be good to have a function which
// periodically checks if some code is no longer minimal, but that would
// need making a global vector with valid DFS code candidates.
// TODO figure out the relationhip between history and rmpath. Why can
// we use the latter to index the former?
void get_min_code_rec(Graph& g, DFSCode& min_dfs, Proj& projected,
                      std::vector<uint32_t>& rmpath) {
  build_rm_path(min_dfs, rmpath);
  // the "from" label of the first edge of the dfs code
  uint32_t min_label = min_dfs[0].from_label;
  // the "to" id of the last forward edge of the dfs code
  uint32_t max_id = min_dfs[rmpath[0]].to;
  Pmap1 root;
  bool added = false; // added new edge?
  uint32_t new_to = 0;
  // generate backward edges
  for (int32_t i = rmpath.size() - 1; !added && i >= 1; --i) {
    // for all edges with the same label triplets that we consider
    // in this function call
    for (int32_t j = 0; j < projected.size(); ++j) {
      PDFS* cur_edge = &projected[j];
      History history;
      build_history(g, cur_edge, history);
      // e1 is expected to be earlier in the DFSCode than e2,
      // because rmpath is built "bottom up".
      Edge *e = get_backward(g, history[rmpath[i]],
                             history[rmpath[0]], history);
      if (e) {
        root[e->elabel].push_back(PDFS(0, e, cur_edge));
        new_to = min_dfs[rmpath[i]].from;
        added = true;
      }
    }
  }
  if (added) {
    min_dfs.push_back(DFSRow(max_id, new_to, g[max_id].label,
                             root.begin()->first, g[new_to].label));
    // the original code has an extra comparison here that I don't need
    get_min_code_rec(g, min_dfs, root.begin()->second, rmpath);
    return;
  }
  Pmap2 root2;
  EdgeList edges;
  uint32_t new_from = 0;
  added = false;
  //generate edges from the last added edge
  for (int32_t j = 0; j < projected.size(); ++j) {
    PDFS* cur_edge = &projected[j];
    History history;
    build_history(g, cur_edge, history);
    if (get_forward_pure(g, history[rmpath[0]],
                         min_label, history, edges)) {
      added = true;
      new_from = max_id;
      for (EdgeList::iterator e = edges.begin();
           e != edges.end(); ++e) {
        root2[(*e)->elabel][g[(*e)->to].label].push_back(
            PDFS(0, *e, cur_edge));
      }
    }
  }
  // generate forward edges
  for (int32_t i = 0; !added && i < rmpath.size(); ++i) {
    for (int32_t j = 0; j < projected.size(); ++j) {
      PDFS* cur_edge = &projected[j];
      History history;
      build_history(g, cur_edge, history);
      if (get_forward_rmpath(g, history[rmpath[i]], min_label,
                             history, edges)) {
        added = true;
        new_from = min_dfs[rmpath[i]].from;
        for (EdgeList::iterator e = edges.begin();
             e != edges.end(); ++e) {
          root2[(*e)->elabel][g[(*e)->to].label].push_back(
              PDFS(0, *e, cur_edge));
        }
      }
    }
  }
  if (added) {
    min_dfs.push_back(DFSRow(new_from, max_id + 1, g[new_from].label,
                             root2.begin()->first,
                             root2.begin()->second.begin()->first));
    // the original code has an extra comparison here that I don't need
    Proj& projected = root2.begin()->second.begin()->second;
    get_min_code_rec(g, min_dfs, projected, rmpath);
    return;
  }
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
void get_min_code(Graph& g, DFSCode& min_dfs) {
  EdgeList edges;
  Pmap3 root;
  std::vector<uint32_t> rmpath;
  min_dfs.clear();
  for (int32_t i = 0; i < g.size(); ++i) {
    get_forward_root(g, g[i], edges);
    for (EdgeList::iterator e = edges.begin(); e != edges.end(); ++e) {
      root[g[i].label][(*e)->elabel][g[(*e)->to].label].push_back(
          PDFS(0, *e, nullptr));
    }
  }
  min_dfs.push_back(
      DFSRow(0, 1, root.begin()->first,
             root.begin()->second.begin()->first,
             root.begin()->second.begin()->second.begin()->first));
  Proj& projected =
      root.begin()->second.begin()->second.begin()->second;
  get_min_code_rec(g, min_dfs, projected, rmpath);
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

void get_min_code_rec(Graph& g, DFSCode& min_dfs, PDFS& cur_edge,
                      std::vector<uint32_t>& rmpath,
                      std::vector<DFSCode>& res_dfs) {
  build_rm_path(min_dfs, rmpath);
  //~ assert(rmpath[0] < min_dfs.size());
  // the "from" label of the first edge of the dfs code
  uint32_t min_label = min_dfs[0].from_label;
  // the "to" id of the last forward edge of the dfs code
  uint32_t max_id = min_dfs[rmpath[0]].to;
  Pmap1 root;
  //~ std::map<uint32_t, std::vector<uint32_t>> new_tos;
  bool added = false; // added new edge?
  uint32_t new_to = 0;
  // generate backward edges
  History history;
  build_history(g, &cur_edge, history);
  for (int32_t i = rmpath.size() - 1; !added && i >= 1; --i) {
    // e1 is expected to be earlier in the DFSCode than e2,
    // because rmpath is built "bottom up".
    //~ assert(rmpath[i] < history.size());
    //~ assert(rmpath[0] < history.size());
    Edge *e = get_backward(g, history[rmpath[i]],
                           history[rmpath[0]], history);
    if (e) {
      root[e->elabel].push_back(PDFS(0, e, &cur_edge));
      //~ new_tos[e->elabel].push_back(min_dfs[rmpath[i]].from);
      new_to = min_dfs[rmpath[i]].from;
      added = true;
    }
  }
  if (added) {
    Proj& projected = root.begin()->second;
    //~ std::vector<uint32_t>& nt = new_tos.begin()->second;
    for (int32_t i = 1; i < projected.size(); ++i) {
      DFSRow new_row = DFSRow(max_id, new_to, g[projected[i].edge->from].label,
                              root.begin()->first, g[projected[i].edge->to].label);
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
      get_min_code_rec(g, min_dfs_local, projected[i], rmpath, res_dfs);
    }
    DFSRow new_row = DFSRow(max_id, new_to, g[projected[0].edge->from].label,
                               root.begin()->first, g[projected[0].edge->to].label);
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
    get_min_code_rec(g, min_dfs, projected[0], rmpath, res_dfs);
    return;
  }
  Pmap2 root2;
  //~ std::map<uint32_t, std::map<uint32_t, std::vector<uint32_t>>> new_froms;
  EdgeList edges;
  uint32_t new_from = 0;
  added = false;
  //generate edges from the last added edge
  //~ History history;
  //~ build_history(g, &cur_edge, history);
    //~ assert(rmpath[0] < history.size());
  if (get_forward_pure(g, history[rmpath[0]],
                       min_label, history, edges)) {
    added = true;
    new_from = max_id;
    for (EdgeList::iterator e = edges.begin();
         e != edges.end(); ++e) {
      root2[(*e)->elabel][g[(*e)->to].label].push_back(
          PDFS(0, *e, &cur_edge));
      new_from = max_id;
      //~ new_froms[(*e)->elabel][g[(*e)->to].label].push_back(max_id);
    }
  }
  // generate forward edges
  for (int32_t i = 0; !added && i < rmpath.size(); ++i) {
    //~ History history;
    //~ build_history(g, &cur_edge, history);
    //~ assert(rmpath[i] < history.size());
    if (get_forward_rmpath(g, history[rmpath[i]], min_label,
                           history, edges)) {
      added = true;
      new_from = min_dfs[rmpath[i]].from;
      for (EdgeList::iterator e = edges.begin();
           e != edges.end(); ++e) {
        root2[(*e)->elabel][g[(*e)->to].label].push_back(
            PDFS(0, *e, &cur_edge));
        //~ new_froms[(*e)->elabel][g[(*e)->to].label].push_back(
            //~ min_dfs[rmpath[i]].from);
        new_from = min_dfs[rmpath[i]].from;
      }
    }
  }
  if (added) {
    Proj& projected = root2.begin()->second.begin()->second;
    //~ std::vector<uint32_t>& nf = new_froms.begin()->second.begin()->second;
    for (int32_t i = 1; i < projected.size(); ++i) {
      //~ assert(i < nf.size());
      DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[i].edge->from].label,
                               root2.begin()->first,
                               root2.begin()->second.begin()->first);
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
      get_min_code_rec(g, min_dfs_local, projected[i], rmpath, res_dfs);
    }
    DFSRow new_row = DFSRow(new_from, max_id + 1, g[projected[0].edge->from].label,
                             root2.begin()->first,
                             root2.begin()->second.begin()->first);
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
    get_min_code_rec(g, min_dfs, projected[0], rmpath, res_dfs);
    return;
  }
  //~ std::cout << "Code candidate:\n";
  //~ print_dfs_code(min_dfs);
  res_dfs.push_back(min_dfs);
}


void get_min_code(Graph& g, std::vector<DFSCode>& res_dfs) {
  EdgeList edges;
  Pmap3 root;
  std::vector<uint32_t> rmpath;
  for (int32_t i = 0; i < g.size(); ++i) {
    get_forward_root(g, g[i], edges);
    for (EdgeList::iterator e = edges.begin(); e != edges.end(); ++e) {
      root[g[i].label][(*e)->elabel][g[(*e)->to].label].push_back(
          PDFS(0, *e, nullptr));
    }
  }
  Proj& projected =
      root.begin()->second.begin()->second.begin()->second;
  for (int32_t i = 0; i < projected.size(); ++i) {
    DFSCode min_dfs_local;
    min_dfs_local.push_back(
      DFSRow(0, 1, root.begin()->first,
             root.begin()->second.begin()->first,
             root.begin()->second.begin()->second.begin()->first));
    get_min_code_rec(g, min_dfs_local, projected[i], rmpath, res_dfs);
  }
}

void generate(Graph& g, Proj& projected, DFSCode& code) {
  DFSCode mincode;
  // TODO check if code is minimal
  
}
