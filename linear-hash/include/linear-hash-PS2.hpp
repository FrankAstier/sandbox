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
        not_split(),
        m(n_buckets),
        b(bucket_size),
        i(0),
        n_entries(0),
        n_probes(0),
        debug(debug)
  {
    for (size_t ii = 0; ii < buckets.size(); ++ii)
      not_split.insert(ii);
  }

  void put(const K& key, const V&value) {

    ++n_entries;
    size_t idx = bucket_index(key, i);
    buckets[idx].push_back({key, value});

    if (buckets[idx].size() > b && not_split.count(idx)) {
      priority.push_back(idx);
      sort(priority.begin(), priority.end(),
           [&](size_t a, size_t b) { return buckets[a].size() < buckets[b].size(); });
    }

    float lf = (float) n_entries / (float) buckets.size();

    if (lf > .95) {

      size_t to_split = numeric_limits<size_t>::max();

      if (priority.back() < (1 << i) * m && not_split.count(idx)) {
        to_split = priority.back();
        priority.resize(priority.size() - 1);
      } else if (!not_split.empty())
        to_split = *not_split.begin();

      assert(to_split != numeric_limits<size_t>::max());

      not_split.erase(to_split);
      size_t new_bucket_idx = to_split + (1 << i) * m;

      if (new_bucket_idx > buckets.size())
        buckets.resize(new_bucket_idx + 1);

      Bucket old_bucket, new_bucket = buckets[new_bucket_idx];

      for (auto kv : buckets[to_split]) {
        size_t h = utils::fasthash64(key);
        size_t k = (1 << (i+1)) * m;
        size_t new_idx = h % k;
        if (new_idx == to_split)
          old_bucket.push_back(kv);
        else
          new_bucket.push_back(kv);
      }

      buckets[idx].swap(old_bucket);

      if (not_split.empty()) {
        ++i;
        assert(buckets.size() == (1 << i) * m);
        for (size_t ii = 0; ii < buckets.size(); ++ii)
          not_split.insert(ii);
      }
    }
  }

  pair<bool, V> get(const K& key) {
    size_t idx = bucket_index(key, i);
    for (size_t j = 0; j < buckets[idx].size(); ++j) {
      if (buckets[idx][j].first == key) {
        n_probes = ((n_entries-1) * n_probes + j)/n_entries;
        return {true, buckets[idx][j].second};
      }
    }
    return {false, V()};
  }

  size_t bucket_index(const K& key, size_t ii) {
    size_t h = utils::fasthash64(key);
    size_t k = (1 << ii) * m;
    size_t hm = h % k;
    if (not_split.count(hm))
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
    << " len:" << n_probes // number of probes to reach a record or not present
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

  vector<Bucket> buckets;
  vector<size_t> priority;
  set<size_t> not_split;
  const size_t m; // initial number of buckets
  const size_t b; // bucket size
  size_t i; // splitting round number
  size_t n_entries; // total number of entries in hash
  float n_probes;
  bool debug;
};

#endif //LINEAR_HASH_PS2_HPP
