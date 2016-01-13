# Perfect hash [![Build status](https://travis-ci.org/FrankAstier/perfect-hash.svg?branch=master)](https://travis-ci.org/FrankAstier/perfect-hash)

## Greedy algorithm
- The first version is a very simple greedy algorithm to break ties (collision in the same bucket) by re-hashing with
  a "salt" that's found greedily, by simply scanning the range of possible salts exhaustively, till we find one
  that simultaneously resolves *all* the collisions in the bucket. That is, we end up with a single salt that works
  to re-hash *all* the keys that were collided in a bucket, to non-colliding buckets. We do that iteratively till
  all the ties are broken.

## References
- Hash, displace and compress: http://cmph.sourceforge.net/papers/esa09.pdf