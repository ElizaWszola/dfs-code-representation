# DFS Code Representation

For now, it's just a small piece of code that supposedly computes the minimum DFS Code of a graph.

Based on GSpan: https://github.com/rkwitt/gboost/tree/master/src-gspan


# Build

```bash
pip install . --user
```

# Run test

```bash
python tests/test_min_dfs.py 
```
# Usage

```python
import _dfs_codes as dfs
# for a graph of n nodes
node_types = [type_1, ..., type_n] # note_types[i] = 'type of vertex with index i'
edge_list = [[from_vertex_id, to_vertex_id, edge_type, edge_id], ...] 
code = dfs.compute_minimal_dfs_code(edge_list, node_types) 
# entries of code look like this: [from_dfs_pos, to_dfs_pos, from_label, edge_label, to_label, from_vertex_id, edge_id, to_vertex_id]
```
