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
 * look at less than 64 bytes, which is a cache line size, and contiguous
 * use fixed size slots for speed: fixed size hash + offset, but has an indirection...
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
    load_factor(0),
    n_probes(0),
    debug(debug)
  {
    for (size_t ii = 0; ii < buckets.size(); ++ii)
      not_split.insert(ii);
  }

  void put(const K& key, const V&value) {

    cout << "Putting: " << key << "," << value << endl;

    ++n_entries;
    size_t idx = bucket_index(key, i);
    if (buckets[idx].empty())
      load_factor += 1.0f/buckets.size();
    buckets[idx].push_back({key, value});

    if (buckets[idx].size() > b && not_split.count(idx) > 0) {
      priority.push_back(idx);
      sort(priority.begin(), priority.end(),
           [&](size_t a, size_t b) { return buckets[a].size() < buckets[b].size(); });
    }

    print_short();

    if (load_factor >= 0.8 || (!priority.empty() && buckets[priority.back()].size() > b)) {

      size_t to_split = numeric_limits<size_t>::max();
      size_t candidate = priority.empty() ? *not_split.begin() : priority.back();

      if (candidate < (1 << i) * m && not_split.count(candidate)) {
        to_split = candidate;
        if (!priority.empty() && to_split == priority.back())
          priority.resize(priority.size() - 1);
      }

      cout << "Splitting: " << to_split << endl;

      assert(to_split != numeric_limits<size_t>::max());

      not_split.erase(to_split);
      size_t new_bucket_idx = to_split + (1 << i) * m;
      buckets.resize(new_bucket_idx + 1);
      Bucket old_bucket;

      cout << "New bucket: " << new_bucket_idx << endl;

      for (auto kv : buckets[to_split]) {
        size_t h = utils::fasthash64(key);
        size_t k = (1 << (i+1)) * m;
        size_t new_idx = h % k;
        if (new_idx == to_split)
          old_bucket.push_back(kv);
        else
          buckets[new_bucket_idx].push_back(kv);
      }

      buckets[to_split].swap(old_bucket);
      if (buckets[to_split].empty())
        load_factor -= 1.0f/buckets.size();

      if (not_split.empty()) {
        ++i;
        for (size_t ii = 0; ii < buckets.size(); ++ii)
          not_split.insert(ii);
      }
    }
  }

  pair<bool, V> get(const K& key) {
    cout << "Retrieving: " << key;
    size_t idx = bucket_index(key, i);
    cout << " in bucket: " << idx << endl;
    print_short();
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
    //cout << " bucket index: " << hm;
    if (not_split.count(hm)) {
      //cout << " bucket has not split, so: " << hm << endl;
      return hm;
    } else {
      //cout << " bucket HAS split, so: " << (h % (2*k)) << endl;
      return h % (2 * k);
    }
  }

  void print_short() const {
    cout << "n buckets: " << buckets.size()
    << " i= " << i
    << " lf= " << load_factor
    << " not split: " << not_split
    << " priority: " << priority
    << endl;
    for (size_t j = 0; j < buckets.size(); ++j) {
      if (debug) cout << "\tBucket " << j;
      if (!buckets[j].empty()) {
        if (debug) cout << " size= " << buckets[j].size() << ": ";
        for (auto kv : buckets[j])
          if (debug) cout << kv << " ";
        if (debug) cout << endl;
      }
      else {
        if (debug) cout << " empty" << endl;
      }
    }
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
  float load_factor;
  float n_probes;
  bool debug;
};

#endif //LINEAR_HASH_PS2_HPP
