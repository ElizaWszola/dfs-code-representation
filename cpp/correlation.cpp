#include "correlation.hpp"

// todo: make a graph out of the code, compute the graph's code, check if equal
bool is_code_min() {
  
  
  return true;
}

void cmc_part_rec(std::vector<Graph>& x, c_real* y, bool* b,
    c_real& best_corr_global, c_real* best_corr_local, uint64_t* local_b,
    double delta, int64_t* max_queue, int64_t* total_processed) {
      

  int64_t local_max_queue = 1;
  int64_t sub_max_queue = 1;
  int64_t local_total_processed = 1;
  
  DFSCode code;
  std::vector<PDFS> candidates;
  std::vector<uint32_t> rmpath;
  std::vector<History> histories;
  for (uint32_t i = 0; i < x.size(); ++i) {
    histories.push_back(History(x[i].size(), x[i].edge_size));
  }
  EdgeList edges;

  c_real mu = std::numeric_limits<c_real>::max();
  cmc_part_rec_subcall(x, y, b, code, candidates, rmpath, histories,
                       edges, best_corr_global,
                       best_corr_local, local_b, delta,
                       mu, local_max_queue, sub_max_queue + 1,
                       local_total_processed);
      
  #if GET_STATS
  (*max_queue) = std::max((*max_queue), local_max_queue);
  (*total_processed) += local_total_processed;
  #endif

}

void cmc_part_rec_subcall(std::vector<Graph>& x, c_real* y, bool* b,
    DFSCode& code, std::vector<PDFS>& candidates,
    std::vector<uint32_t>& rmpath, std::vector<History>& histories,
    EdgeList& edges, c_real& best_corr_global, c_real* best_corr_local,
    uint64_t* local_b, double delta, c_real curr_mu,
    int64_t& max_queue, int64_t sub_max_queue, int64_t& total_processed) {
      
  c_real mu1 = 0.0;
  c_real mu2 = 0.0;
  
  #if GET_STATS
  max_queue = std::max(max_queue, sub_max_queue);
  ++total_processed;
  #endif
  
  c_real bcl = *best_corr_local;
  c_real bcg = best_corr_global;
  if (curr_mu <= (std::fabs(bcg) > std::fabs(bcl) ? std::fabs(bcg) : std::fabs(bcl)) + delta) {
    return;
  }
 
  // add to mu1/mu2 for each graph the subgraph is a part of
  uint32_t last_gid = std::numeric_limits<uint32_t>::max();
  for (int32_t i = 0; i < candidates.size(); ++i) {
    uint32_t pid = candidates[i].id;
    if (pid != last_gid) {
      last_gid = pid;
      c_real y1 = y[pid];
      bool ycz = (y1 > 0);
      mu1 += ycz * y1;
      mu2 -= !ycz * y1;
    }
  }
  
  c_real corr = mu1 - mu2;
  bcg = best_corr_global;
  c_real bc = std::fabs(bcg) > std::fabs(bcl) ? std::fabs(bcg) : std::fabs(bcl);
  if (std::fabs(corr) > bc) {
    *best_corr_local = corr;
    bc = std::fabs(corr);
    // TODO: copy the best solution
    //~ std::memcpy(local_b, v, featbits * sizeof(uint64_t));
  }
  c_real max_mu = std::max(mu1, mu2);
  if (max_mu > bc + delta) {
  
    uint32_t rmp_size = build_rm_path(code, rmpath);
    uint32_t min_label = code[0].from_label;
    uint32_t max_id = code[rmpath[0]].to;
    Pmap3 fwd_root;
    Pmap2 bwd_root;

    // generate new subgraphs into maps
    
    for (uint32_t cid = 0; cid < candidates.size(); ++cid) {
      
      PDFS& cur_edge = candidates[cid];
      uint32_t gid = cur_edge.id;
      Graph& g = x[gid];
      History& history = histories[gid];
      
      uint32_t hsize = build_history(g, &cur_edge, history);
      
      for (int32_t i = rmp_size - 1; i >= 1; --i) {
        Edge *e = get_backward(g, history[rmpath[i]],
                               history[rmpath[0]], history);
        if (e) {
          bwd_root[g[e->from].label][e->elabel].push_back(PDFS(gid, e, &cur_edge));
        }
      }
      //generate edges from the last added edge
      uint32_t edge_size = get_forward_pure(g, history[rmpath[0]],
                                            min_label, history, edges);
      for (uint32_t eid = 0; eid < edge_size; ++eid) {
        Edge* e = edges[eid];
        fwd_root[g[e->from].label][e->elabel][g[e->to].label].push_back(
            PDFS(gid, e, &cur_edge));

        for (int32_t i = 0; i < rmp_size; ++i) {
          uint32_t edge_size = get_forward_rmpath(g, history[rmpath[i]],
                                                  min_label, history, edges);
          for (uint32_t eid = 0; eid < edge_size; ++eid) {
            Edge* e = edges[eid];
            fwd_root[g[e->from].label][e->elabel][g[e->to].label].push_back(
                PDFS(gid, e, &cur_edge));
          }
        }
      }
    }
        
      // loop through the maps
    for (Pmap2::iterator to = bwd_root.begin(); to != bwd_root.end(); ++to) {
      for (Pmap1::iterator elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
        code.push_back(DFSRow(max_id, to->first, -1, elabel->first, -1));
        cmc_part_rec_subcall(x, y, b, code, elabel->second, rmpath,
                             histories, edges, best_corr_global, best_corr_local,
                             local_b, delta, curr_mu, max_queue,
                             sub_max_queue, total_processed);
        code.pop_back();
      }
    }
    for (Pmap3::iterator from = fwd_root.begin(); from != fwd_root.end(); ++from) {
      for (Pmap2::iterator elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
        for (Pmap1::iterator to = elabel->second.begin(); to != elabel->second.end(); ++to) {
          code.push_back(DFSRow(from->first, max_id + 1, -1, elabel->first, -1));
          cmc_part_rec_subcall(x, y, b, code, to->second, rmpath,
                             histories, edges, best_corr_global, best_corr_local,
                             local_b, delta, curr_mu, max_queue,
                             sub_max_queue, total_processed);
          code.pop_back();
        }
      }
    }
  }
}
