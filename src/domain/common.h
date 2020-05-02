#ifndef COMMON
#define COMMON
#include <cstddef>
#include <unordered_map>
#include <utility>
using namespace std;
struct pair_hash {
  template <typename T1, typename T2> size_t operator()(const pair<T1, T2> &pair) const {
    return hash<T1>()(pair.first) ^ hash<T2>()(pair.second);
  }
};

typedef unordered_map<pair<string, string>, double, pair_hash> LabelPairMap;
#endif /* COMMON */
