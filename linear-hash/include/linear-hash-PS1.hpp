/**
 * Linear hashing prototype.
 */

#ifndef LINEAR_HASH_PS1_HPP
#define LINEAR_HASH_PS1_HPP

#include <fast_hash.hpp>
#include <stl_io.hpp>

using namespace utils;

/**
 * Linear hash with priority splitting - take 1.
 *
 * Here, we opportunistically split a bucket that would overflow with the current
 * insertion, if there is already a bucket created to accommodate some of the keys
 * (hi and hi+1 send to bucket j or bucket j + m always).
 */
template <typename K, typename V>
struct LinearHash_PS1 {

  LinearHash_PS1(size_t n_buckets =100, size_t bucket_size =4, bool debug = false)
      : buckets(n_buckets),
        priority(),
        m(n_buckets),
        b(bucket_size),
        i(0),
        p(0),
        n_entries(0),
        n_hits(0),
        search_length(0),
        debug(debug)
  {}

  void put(const K& key, const V&value) {
    if (debug) std::cout << "Inserting " << key << " " << value << " - ";
    if (debug) std::cout << "m= " << m << " b= " << b << " i= " << i << " p= " << p;
    if (debug) std::cout << " n buckets= " << buckets.size() << std::endl;
    ++n_entries;
    uint64_t idx = bucket_index(key, i);
    assert(idx < buckets.size());
    assert(p <= (1 << i) * m);
    assert(buckets.size() == (1 << i) * m + p);
    if (debug) std::cout << "\tStoring " << key << "," << value << " in bucket: " << idx << std::endl;
    buckets[idx].push_back({key, value});

    if (buckets[idx].size() > b && idx < p) {
      ++n_hits;
      if (debug) std::cout << "\tOverflow on bucket " << idx << " - Splitting bucket " << idx << std::endl;
      Bucket old_bucket;
      for (auto kv : buckets[idx]) {
        uint64_t h = utils::fasthash64(key);
        uint64_t k = (1 << (i+1)) * m;
        uint64_t new_idx = h % k;
        Bucket new_bucket = buckets[new_idx];
        if (debug) std::cout << "\tkey: " << kv.first;
        if (debug) std::cout << " old index= " << bucket_index(kv.first, i);
        if (debug) std::cout << " -> new index= " << new_idx << " ";
        assert(new_idx <= buckets.size());
        assert(new_idx == idx || new_idx == idx + (1 << i) * m);
        if (new_idx == idx) {
          if (debug) std::cout << "\t\t" << kv << " stays in bucket " << idx << std::endl;
          old_bucket.push_back(kv);
        } else {
          if (debug) std::cout << "\t\t" << kv << " goes to new bucket " << (idx+m) << std::endl;
          new_bucket.push_back(kv);
        }
      }
      buckets[idx].swap(old_bucket);

    } else if (buckets[idx].size() > b) {
      if (debug) std::cout << "\tOverflow on bucket " << idx << " - Splitting bucket " << p << std::endl;
      Bucket old_bucket, new_bucket;
      for (auto kv : buckets[p]) {
        uint64_t new_idx = bucket_index(kv.first, i+1);
        if (debug) std::cout << "\tnew index= " << new_idx << " ";
        assert(new_idx <= buckets.size());
        assert(new_idx == p || new_idx == p + (1 << i) * m);
        if (new_idx == p) {
          if (debug) std::cout << "\t\t" << kv << " stays in bucket " << p << std::endl;
          old_bucket.push_back(kv);
        } else {
          if (debug) std::cout << "\t\t" << kv << " goes to new bucket " << new_idx << std::endl;
          new_bucket.push_back(kv);
        }
      }
      buckets[p].swap(old_bucket);
      buckets.push_back(new_bucket);
      ++p;
      if (p == (1 << i) * m) {
        if (debug) std::cout << "\tp = " << p << " going to next round" << std::endl;
        p = 0;
        ++i;
      }
    }
  }

  std::pair<bool, V> get(const K& key) {
    uint64_t idx = bucket_index(key, i);
    for (size_t j = 0; j < buckets[idx].size(); ++j) {
      if (buckets[idx][j].first == key) {
        search_length = ((n_entries-1) * search_length + j)/n_entries;
        return {true, buckets[idx][j].second};
      }
    }
    return {false, V()};
  }

  uint64_t bucket_index(const K& key, size_t ii) const {
    uint64_t h = utils::fasthash64(key);
    uint64_t k = (1 << ii) * m;
    uint64_t hm = h % k;
    uint64_t h2 = h % (2*k);
    return hm >= p ? hm : h2;
  }

  void print() const {
    size_t n_empty_buckets = 0;
    float avg_bucket_occupancy = 0;
    float min_bucket_occupancy = std::numeric_limits<float>::max();
    float max_bucket_occupancy = 0;

    for (size_t j = 0; j < buckets.size(); ++j) {
      if (debug) std::cout << "\tBucket " << j;
      if (!buckets[j].empty()) {
        if (debug) std::cout << " size= " << buckets[j].size() << ": ";
        for (auto kv : buckets[j])
          if (debug) std::cout << kv << " ";
        if (debug) std::cout << std::endl;
        min_bucket_occupancy = std::min(min_bucket_occupancy, (float) buckets[j].size());
        max_bucket_occupancy = std::max(max_bucket_occupancy, (float) buckets[j].size());
        avg_bucket_occupancy += buckets[j].size();
      }
      else {
        if (debug) std::cout << " empty" << std::endl;
        ++n_empty_buckets;
      }
    }

    std::cout << "n buckets:" << buckets.size()
    << " entries:" << n_entries
    << " empty:" << n_empty_buckets
    << " min:" << min_bucket_occupancy
    << " avg:" << (avg_bucket_occupancy/(buckets.size() - n_empty_buckets))
    << " max:" << max_bucket_occupancy
    << " hits:" << n_hits
    << " len:" << search_length
    << std::endl;

    std::vector<float> dist((size_t)max_bucket_occupancy+1);

    for (Bucket bucket : buckets)
      ++dist[bucket.size()];
    float M = 0;
    for (size_t j = 0; j < dist.size(); ++j) {
      dist[j] /= (float) n_entries;
      M = std::max(M, dist[j]);
    }
    std::cout << "\ndist= " << dist << std::endl << std::endl;
    float k = 30 / M;
    for (size_t j = 0; j < dist.size(); ++j) {
      size_t l =  (size_t) (k * dist[j]);
      std::cout << j << "\t:";
      for (size_t jj = 0; jj < l; ++jj)
        std::cout << "*";
      std::cout << std::endl;
    }
  }

  typedef std::vector<std::pair<K,V>> Bucket;

  std::vector<Bucket> buckets;
  std::vector<size_t> priority;
  size_t m; // initial number of buckets
  size_t b; // bucket size
  size_t i; // splitting round number
  size_t p; // bucket to split next
  size_t n_entries; // total number of entries in hash
  size_t n_hits; // lucky hits where we split the right bucket
  float search_length;
  bool debug;
};

#endif //LINEAR_HASH_PS1_HPP
