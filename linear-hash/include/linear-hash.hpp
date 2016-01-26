/**
 * Linear hashing prototype.
 */

#ifndef LINEAR_HASH_HPP
#define LINEAR_HASH_HPP

#include <fast_hash.hpp>

template <typename K, typename V>
struct LinearHash {

  LinearHash(size_t n_buckets =100)
  : buckets(n_buckets),
    i(0),
    p(0)
  {}

  void put(const K& key, const V& bytes) {

  }

  V get(const std::string& key) const {
    return V();
  }

  uint64_t bucket_index(const K& key) const {
    uint64_t h = utils::fasthash64(key);
    uint64_t k = (1 << i) * buckets.size();
    uint64_t hm = h % k;
    return hm >= p ? hm : h % (2*k);
  }

  typedef std::vector<std::pair<K,V>> Bucket;

  std::vector<Bucket> buckets;
  size_t i; // splitting round number
  size_t p; // bucket to split next
};

#endif //LINEAR_HASH_HPP
