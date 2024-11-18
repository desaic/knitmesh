#ifndef EDGE_H
#define EDGE_H
struct Edge {
  unsigned v[2];
  // sorts the indice increasing order
  Edge(unsigned v0, unsigned v1) {
    if (v0 <= v1) {
      v[0] = v0;
      v[1] = v1;
    } else {
      v[0] = v1;
      v[1] = v0;
    }
  }
  bool operator<(const Edge& o) const {
    if (v[0] < o.v[0]) {
      return true;
    }
    return (v[0] == o.v[0] && v[1] < o.v[1]);
  }

  bool operator==(const Edge& o) const {
    return v[0] == o.v[0] && v[1] == o.v[1];
  }
};
#endif