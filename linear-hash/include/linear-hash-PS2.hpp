/**
 * Linear hashing prototype.
 */

#ifndef LINEAR_HASH_PS2_HPP
#define LINEAR_HASH_PS2_HPP

#include <fast_hash.hpp>
#include <stl_io.hpp>

using namespace utils;
using namespace std;

/**
 * Linear hash with priority splitting - take 2: splitting with priority.
 *
 * linear probe
 * split or fill holes in zone 3? just split first
 */
template <typename K, typename V>
struct LinearHash_PS2 {

  LinearHash_PS2(size_t n_buckets =100, size_t bucket_size =4, bool debug = false)
  : buckets(n_buckets),
    priority(),
    not_split(buckets.size(), 1),
    m(n_buckets),
    b(bucket_size),
    i(0),
    n_entries(0),
    search_length(0),
    debug(debug)
  {}

  void put(const K& key, const V&value) {

//    ++n_entries;
//    size_t idx = bucket_index(key, i);
//
//    buckets[idx].push_back({key, value});
//
//    if (buckets[idx].size() > b && not_split[idx])
//      priority.push_back(idx);
//
//    sort(priority.begin(), priority.end(), order);
//
//    float lf = (float) n_entries / (float) buckets.size();
//
//    if (lf > .95) {
//
//      size_t to_split = priority.back() < (1 << i) * m ? priority.back() : *not_split.begin();
//      priority.resize(priority.size() - 1);
//      not_split.erase(to_split);
//
//      size_t new_idx = to_split + (1 << i) * m;
//
//      if (new_idx > buckets.size())
//        buckets.resize(new_idx + 1);
//
//      Bucket old_bucket, new_bucket = buckets[new_idx];
//
//      for (auto kv : buckets[to_split]) {
//
//        if (new_idx == to_split)
//          old_bucket.push_back(kv);
//        else
//          new_bucket.push_back(kv);
//      }
//
//      buckets[idx].swap(old_bucket);
//
//      p = *not_split.begin();
//      if (p == (1 << i) * m) {
//        p = 0;
//        ++i;
//        size_t new_size = (1 << (i+1)) * m;
//        not_split.resize(new_size);
//        fill(not_split.begin(), not_split.end(), 1);
//      }
//    }
  }

  pair<bool, V> get(const K& key) {
    size_t idx = bucket_index(key, i);
    for (size_t j = 0; j < buckets[idx].size(); ++j) {
      if (buckets[idx][j].first == key) {
        search_length = ((n_entries-1) * search_length + j)/n_entries;
        return {true, buckets[idx][j].second};
      }
    }
    return {false, V()};
  }

  size_t bucket_index(const K& key, size_t ii) const {
    size_t h = utils::fasthash64(key);
    size_t k = (1 << ii) * m;
    size_t hm = h % k;
    if (not_split[hm])
      return hm;
    return h % (2*k);
  }

  void print() const {
    size_t n_empty_buckets = 0;
    float avg_bucket_occupancy = 0;
    float min_bucket_occupancy = numeric_limits<float>::max();
    float max_bucket_occupancy = 0;

    for (size_t j = 0; j < buckets.size(); ++j) {
      if (debug) cout << "\tBucket " << j;
      if (!buckets[j].empty()) {
        if (debug) cout << " size= " << buckets[j].size() << ": ";
        for (auto kv : buckets[j])
          if (debug) cout << kv << " ";
        if (debug) cout << endl;
        min_bucket_occupancy = min(min_bucket_occupancy, (float) buckets[j].size());
        max_bucket_occupancy = max(max_bucket_occupancy, (float) buckets[j].size());
        avg_bucket_occupancy += buckets[j].size();
      }
      else {
        if (debug) cout << " empty" << endl;
        ++n_empty_buckets;
      }
    }

    cout << "n buckets:" << buckets.size()
    << " entries:" << n_entries
    << " empty:" << n_empty_buckets
    << " min:" << min_bucket_occupancy
    << " avg:" << (avg_bucket_occupancy/(buckets.size() - n_empty_buckets))
    << " max:" << max_bucket_occupancy
    << " len:" << search_length // number of probes to reach a record or not present
    << endl;

    vector<float> dist((size_t)max_bucket_occupancy+1);

    for (Bucket bucket : buckets)
      ++dist[bucket.size()];
    float M = 0;
    for (size_t j = 0; j < dist.size(); ++j) {
      dist[j] /= (float) n_entries;
      M = max(M, dist[j]);
    }
    cout << "\ndist= " << dist << endl << endl;
    float k = 30 / M;
    for (size_t j = 0; j < dist.size(); ++j) {
      size_t l =  (size_t) (k * dist[j]);
      cout << j << "\t:";
      for (size_t jj = 0; jj < l; ++jj)
        cout << "*";
      cout << endl;
    }
  }

  typedef vector<pair<K,V>> Bucket;

  function<bool(size_t,size_t)> order = [&](size_t a, size_t b) {
    return buckets[a].size() < buckets[b].size();
    //return !((a < p && b >= p)
    //       || (a < p && b < p && buckets[a].size() > buckets[b].size()));
  };

  vector<Bucket> buckets;
  vector<size_t> priority;
  vector<bool> not_split;
  const size_t m; // initial number of buckets
  const size_t b; // bucket size
  size_t i; // splitting round number
  size_t n_entries; // total number of entries in hash
  float search_length;
  bool debug;
};

#endif //LINEAR_HASH_PS2_HPP
