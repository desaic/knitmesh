#pragma once
#include <vector>
#include <set>

/// @brief disjoint set union find.
class DisjSet {
 public:
  struct UFNode {
    UFNode(unsigned p = 0) : parent(p) {}
    unsigned parent = 0;
    // self 1 + number of children
    unsigned size = 1;
  };
  /// @brief Create n roots.
  void Create(unsigned numNodes) { 
    nodes.resize(numNodes);
    for (unsigned i = 0; i < nodes.size(); i++) {
      nodes[i].parent = i;
    }
  }
  /// @brief Merge root of i and root of j by size.
  /// @return root of merged tree
  unsigned Union(unsigned i, unsigned j) {
    if (i >= nodes.size() || j >= nodes.size()) {
      return 0;
    }
    unsigned rooti = Find(i);
    unsigned rootj = Find(j);
    if (rooti == rootj) return rooti;
    unsigned sizei = nodes[rooti].size;
    unsigned sizej = nodes[rootj].size;
    if (sizei<sizej) {
      nodes[rooti].parent = rootj;
      nodes[rootj].size += nodes[rooti].size;
      return rootj;
    } else {
      nodes[rootj].parent = rooti;
      nodes[rooti].size += nodes[rootj].size;
      return rooti;
    }
  }
  /// find root node for i
  unsigned Find(unsigned i) {
    if (i >= nodes.size()) {
      return 0;
    }
    if (nodes[i].parent == i) {
      return i;
    } else {
      unsigned root = Find(nodes[i].parent);
      nodes[i].parent = root;
      return root;
    }
  }
  /// find root node for i const version
  unsigned Find(unsigned i) const {
    if (i >= nodes.size()) {
      return 0;
    }
    if (nodes[i].parent == i) {
      return i;
    } else {
      unsigned root = Find(nodes[i].parent);
      return root;
    }
  }
  bool IsRoot(unsigned i) { return i < nodes.size() && nodes[i].parent == i; }
  std::vector<unsigned> GetRoots() const { std::set<unsigned> roots;
    for (unsigned i = 0; i < nodes.size(); i++) {
      roots.insert(Find(i));
    }
    return std::vector<unsigned>(roots.begin(), roots.end());
  }
 private:
  std::vector<UFNode> nodes;
};
