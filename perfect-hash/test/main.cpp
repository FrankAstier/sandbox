#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <set>

#include <fast_hash.hpp>
#include <stl_io.hpp>

using namespace std;
using namespace utils;

using Bucket = list<string>;

vector<string> generate_keys(size_t N) {

  vector<string> keys(N);

  for (size_t i = 0; i < N; ++i) {
    stringstream s;
    s << "key-" << i;
    keys[i] = s.str();
  }

  return keys;
}

/**
 * Greedy algorithm to find salts for perfect hashing.
 *
 * Notes:
 * - could leave some buckets with > 1 items
 * - minimum salt size?
 */
vector<int> find_H(const vector<string>& keys, bool _debug_algo =false) {

  size_t N = keys.size();
  vector<int> H(N);
  vector<Bucket> buckets(N);

  // Hash to buckets
  for (size_t i = 0; i < N; ++i) {
    size_t h0 = fasthash64(keys[i]) % N;
    buckets[h0].push_back(keys[i]);
  }

  // Find salts
  sort(buckets.begin(), buckets.end(),
       [](const Bucket& b1, const Bucket& b2) { return b1.size() > b2.size(); });

  size_t M = numeric_limits<size_t>::max();
  size_t last = 0;
  set<size_t> free;
  size_t n = 0;
  generate_n(inserter(free, free.end()), N, [&]() { return n++;});

  if (_debug_algo) cout << buckets << endl;

  for (Bucket bucket : buckets) {
    if (bucket.size() <= 1) break;
    ++last;
    set<size_t> slots;
    size_t salt = 1;
    size_t h0 = fasthash64(bucket.front()) % N;
    if (_debug_algo) cout << "Bucket: " << bucket << " @ " << h0 << endl;
    for (; salt < M && slots.size() < bucket.size(); ++salt) {
      if (_debug_algo) cout << "\tTrying salt: " << salt << endl;
      if (_debug_algo) cout << "\t\tfree: " << free << " slots: " << slots << endl;
      for (auto key : bucket) {
        size_t h1 = fasthash64(key, salt) % N;
        if (free.find(h1) == free.end()) {
          if (_debug_algo) cout << "\t\tkey: " << key << " slot: " << h1 << " - already occupied" << endl;
          free.insert(slots.begin(), slots.end()); // backtrack
          slots.clear();
          break;
        } else {
          slots.insert(h1);
          free.erase(h1);
          if (_debug_algo) cout << "\t\tkey: " << key << " taking slot: " << h1 << endl;
          if (_debug_algo) cout << "\t\tfree: " << free << " slots: " << slots << endl;
        }
      }
    }
    if (slots.size() == bucket.size()) {
      size_t h0 = fasthash64(bucket.front()) % N;
      H[h0] = (int) salt-1;
      if (_debug_algo) cout << "\tGood salt: " << (salt-1) << " stored at: " << h0 << endl;
    } else {
      if (_debug_algo) cout << "\tInfeasible!" << endl;
    }
  }

  // Handle size 1 buckets if any
  for (size_t i = last; i < buckets.size(); ++i) {
    if (buckets[i].empty()) break;
    if (_debug_algo) cout << "Size 1 bucket: " << buckets[i] << endl;
    size_t h0 = fasthash64(buckets[i].front()) % N;
    int f = (int) *free.begin();
    free.erase(*free.begin());
    H[h0] = -f - 1;
    if (_debug_algo) cout << "\t@ " << f << endl;
    if (_debug_algo) cout << "\tfree: " << free << endl;
  }

  return H;
}

size_t p_hash(const string& key, const vector<int>& H) {
  size_t h0 = fasthash64(key) % H.size();
  return (size_t) (H[h0] < 0 ? (size_t) - H[h0] - 1 : fasthash64(key, (size_t) H[h0])) % H.size();
}

void test_perfect_hash() {

  size_t M = 3, N = 1000;

  for (size_t i = 0; i < M; ++i) {
    vector<string> keys = generate_keys(N);
    vector<int> H = find_H(keys);

    set<size_t> results;
    for (string key : keys) {
      size_t h = p_hash(key, H);
      if (results.count(h)) {
        cout << "Collision seen " << key << " @ " << h << endl;
        exit(-9);
      }
      results.insert(h);
    }
    cout << "Test successful " << i << " - no collision" << endl;
  }
}

int main() {

  test_perfect_hash();
  return 0;
}