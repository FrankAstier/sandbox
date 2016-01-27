/**
 * Linear hashing prototype.
 */

#ifndef LINEAR_HASH_HPP
#define LINEAR_HASH_HPP

#include <fast_hash.hpp>
#include <stl_io.hpp>

using namespace utils;

template <typename K, typename V>
struct LinearHash {

  LinearHash(size_t n_buckets =100, size_t bucket_size =4)
      : buckets(n_buckets),
        m(n_buckets),
        b(bucket_size),
        i(0),
        p(0)
  {}

  void put(const K& key, const V&value) {
    std::cout << "Inserting " << key << " " << value << " - ";
    std::cout << "m= " << m << " b= " << b << " i= " << i << " p= " << p;
    std::cout << " n buckets= " << buckets.size() << std::endl;
    uint64_t idx = bucket_index(key, i);
    assert(idx < buckets.size());
    assert(p <= (1 << i) * m);
    std::cout << "\tStoring " << key << "," << value << " in bucket: " << idx << std::endl;
    buckets[idx].push_back({key, value});
    if (buckets[idx].size() > b) { // overflow
      std::cout << "\tOverflow on bucket " << idx << " - Splitting bucket " << p << std::endl;
      Bucket old_bucket, new_bucket;
      for (auto kv : buckets[p]) {
        uint64_t new_idx = bucket_index(kv.first, i+1);
        std::cout << "\tnew index= " << new_idx << " ";
        assert(new_idx <= buckets.size() + 1);
        assert(new_idx == p || new_idx == p + (1 << i) * m);
        if (new_idx == p) {
          std::cout << "\t\t" << kv << " stays in bucket " << idx << std::endl;
          old_bucket.push_back(kv);
        } else {
          std::cout << "\t\t" << kv << " goes to new bucket " << (idx+m) << std::endl;
          new_bucket.push_back(kv);
        }
      }
      buckets[p].swap(old_bucket);
      buckets.push_back(new_bucket);
      ++p;
      if (p == (1 << i) * m) {
        std::cout << "\tp = " << p << " going to next round" << std::endl;
        p = 0;
        ++i;
      }
    }
  }

  std::pair<bool, V> get(const K& key) const {
    uint64_t idx = bucket_index(key, i);
    for (auto kv : buckets[idx])
      if (kv.first == key)
        return {true, kv.second};
    return {false, V()};
  }

  uint64_t bucket_index(const K& key, size_t ii) const {
    uint64_t h = utils::fasthash64(key);
    uint64_t k = (1 << ii) * m;
    uint64_t hm = h % k;
    return hm >= p ? hm : h % (2*k);
  }

  void print() const {
    for (size_t j = 0; j < buckets.size(); ++j) {
      std::cout << "\tBucket " << j;
      if (!buckets[j].empty()) {
        std::cout << " size= " << buckets[j].size() << ": ";
        for (auto kv : buckets[j])
          std::cout << kv << " ";
        std::cout << std::endl;
      }
      else
        std::cout << " empty" << std::endl;
    }
  }

  typedef std::vector<std::pair<K,V>> Bucket;

  std::vector<Bucket> buckets;
  size_t m; // initial number of buckets
  size_t b; // bucket size
  size_t i; // splitting round number
  size_t p; // bucket to split next
};

#endif //LINEAR_HASH_HPP
