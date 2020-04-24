#include "../domain/graph.h"

template <typename N, typename E> class GXL {

public:
  void write(const Graph<N, E> &g, const string &fname) const;
};