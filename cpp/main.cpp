#ifdef FAST
#include "graph_fast.hpp"
#else
#include "graph.hpp"
#endif
#include <iostream>

int main() {
  std::cout << "Computing minimal DFS Code" << std::endl;
  Graph g;
  DFSCode result;
  create_sample(g);
  get_min_code(g, result, 0);
  //~ print_graph(g);
  std::vector<DFSCode> results;
  //~ get_min_code(g, results);
  print_dfs_code(result);
  for (int32_t i = 0; i < results.size(); ++i) {
    std::cout << "Min DFS Candidate #" << i + 1 << "\n";
    print_dfs_code(results[i]);
    std::cout << "\n";
  }
  return 0;
}
