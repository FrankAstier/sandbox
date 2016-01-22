/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2015 Frank Astier
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef CPP_UTILS_STL_IO_HPP
#define CPP_UTILS_STL_IO_HPP

#include <iostream>
#include <iomanip>
#include <iterator>

#include <type_traits.hpp>

namespace utils {

/**
 * A class to manipulate the way we display STL objects.
 */
struct stl_io_manips {
  static size_t abbreviate_threshold; // 0 means "don't abbreviate"
  static bool abbreviate_display_n; // display number of elements left
};

size_t stl_io_manips::abbreviate_threshold = 0;
bool stl_io_manips::abbreviate_display_n = false;

struct SetAbbreviateThreshold { size_t t; };
inline SetAbbreviateThreshold set_abbreviate_threshold(size_t tt) { return {tt}; }

inline std::ostream& operator<<(std::ostream& out, SetAbbreviateThreshold sat) {
  stl_io_manips::abbreviate_threshold = sat.t;
  return out;
}

struct SetAbbreviateDisplayN { bool display; };
inline SetAbbreviateDisplayN set_abbreviate_display_n(bool display =true) { return {display}; }

inline std::ostream& operator<<(std::ostream& out, SetAbbreviateDisplayN sadn) {
  stl_io_manips::abbreviate_display_n = sadn.display;
  return out;
}

/**
 * Stream a pair. No surprise here.
 */
template <typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& out, const std::pair<T1,T2>& p) {
  return out << "(" << p.first << "," << p.second << ")";
}

/**
 * Dump the contents of any std container to a stream. Very useful in debugging
 * in particular.
 */
template <typename C, typename std::enable_if<utils::is_std_container<C>::value>::type* =nullptr>
inline std::ostream& operator<<(std::ostream& out, const C& c) {

  size_t S = c.size();

  if (S == 0)
    return out << "{}";

  out << "{";

  if (stl_io_manips::abbreviate_threshold > 0)
    S = std::min(S, stl_io_manips::abbreviate_threshold+1);

  auto it = c.begin();
  for (size_t i = 0; i < S - 1; ++i, ++it)
    out << *it << ",";

  if (stl_io_manips::abbreviate_threshold > 0)
    if (stl_io_manips::abbreviate_display_n) {
      out << "..(" << c.size()-(S-1) << ")..";
    } else {
      out << "..";
    }
  else
    out << *it;

  out << "}";
  return out;
}

///**
// * Stream integer in binary, with group and leading zeros if needed.
// */
//inline std::string bin(unsigned long long x, size_t g =8, size_t n =64) {
//  size_t n_bits = (size_t) std::log2(x) + 1;
//  std::cout << n_bits << std::endl;
//  bool str[n];
//  for (size_t i = 0; i < n_bits; ++i) {
//    str[n-i-1] = (bool) ((x & (1 << i)) / (1 << i));
//  }
//  std::stringstream s;
//  for (size_t i = 0; i < n - n_bits; ++i) {
//    s << '0';
//    if (i % g == 0)
//      s << ' ';
//  }
//  for (size_t i = 0; i < n_bits; ++i) {
//    s << str[i];
//    if (i % g == 0)
//      s << ' ';
//  }
//  return s.str();
//}

////--------------------------------------------------------------------------------
//// IO CONTROL AND MANIPULATORS
////--------------------------------------------------------------------------------
//struct IOControl
//{
//  int abbr;                  // shorten long vectors output
//  bool output_n_elts;        // output vector size at beginning
//
//  bool pair_paren;           // put parens around pairs in vector of pairs
//  const char* pair_sep;      // put separator between pair.first and pair.second
//
//  bool bit_vector;           // output 0/1 vector compactly
//
//  inline IOControl(int a =-1, bool s =true, bool pp =false, const char* psep =" ",
//                   bool bv =false)
//  : abbr(a),
//    output_n_elts(s),
//    pair_paren(pp),
//    pair_sep(psep),
//    bit_vector(bv)
//  {}
//
//  inline void reset()
//  {
//    abbr = -1;
//    output_n_elts = true;
//    pair_paren = false;
//    pair_sep = " ";
//    bit_vector = false;
//  }
//};
//
//extern IOControl io_control;
//
//template <typename CharT, typename Traits, typename T>
//inline std::basic_ostream<CharT,Traits>&
//operator,(std::basic_ostream<CharT,Traits>& out_stream, const T& a)
//{
//  return out_stream << ' ' << a;
//}
//
//template <typename CharT, typename Traits, typename T>
//inline std::basic_istream<CharT,Traits>&
//operator,(std::basic_istream<CharT,Traits>& in_stream, T& a)
//{
//  return in_stream >> a;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//operator,(std::basic_ostream<CharT,Traits>& out_stream,
//          std::basic_ostream<CharT,Traits>& (*pf)(std::basic_ostream<CharT,Traits>&))
//{
//  pf(out_stream);
//  return out_stream;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//p_paren(std::basic_ostream<CharT,Traits>& out_stream)
//{
//  io_control.pair_paren = true;
//  return out_stream;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//psep_comma(std::basic_ostream<CharT,Traits>& out_stream)
//{
//  io_control.pair_sep = ",";
//  return out_stream;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//psep_dot(std::basic_ostream<CharT,Traits>& out_stream)
//{
//  io_control.pair_sep = ".";
//  return out_stream;
//}
//
//struct abbr
//{
//  int n;
//  inline abbr(int _n) : n(_n) {}
//};
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//operator<<(std::basic_ostream<CharT,Traits>& out_stream, abbr s)
//{
//  io_control.abbr = s.n;
//  return out_stream;
//}
//
//struct debug
//{
//  int n;
//  inline debug(int _n =-1) : n(_n) {}
//};
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//operator<<(std::basic_ostream<CharT,Traits>& out_stream, debug d)
//{
//  io_control.abbr = d.n;
//  io_control.output_n_elts = false;
//  io_control.pair_sep = ",";
//  io_control.pair_paren = true;
//  return out_stream;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_istream<CharT,Traits>&
//bit_vector(std::basic_istream<CharT,Traits>& in_stream)
//{
//  io_control.bit_vector = true;
//  return in_stream;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//bit_vector(std::basic_ostream<CharT,Traits>& out_stream)
//{
//  io_control.bit_vector = true;
//  return out_stream;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_istream<CharT,Traits>&
//general_vector(std::basic_istream<CharT,Traits>& in_stream)
//{
//  io_control.bit_vector = false;
//  return in_stream;
//}
//
//template <typename CharT, typename Traits>
//inline std::basic_ostream<CharT,Traits>&
//general_vector(std::basic_ostream<CharT,Traits>& out_stream)
//{
//  io_control.bit_vector = false;
//  return out_stream;
//}
//
////--------------------------------------------------------------------------------
//// CHECKERS
////--------------------------------------------------------------------------------
//template <typename T1>
//struct is_positive_checker
//{
//  T1& var;
//
//  inline is_positive_checker(T1& v) : var(v) {}
//
//  template <typename CharT, typename Traits>
//  inline void do_check(std::basic_istream<CharT,Traits>& in_stream)
//  {
//    double value = 0;
//    in_stream >> value;
//    if (value < 0) {
//      std::cout << "Value out of range: " << value
//      << " - Expected positive or zero value"
//      << std::endl;
//      exit(-1);
//    }
//    var = (T1) value;
//  }
//};
//
//template <typename CharT, typename Traits, typename T1>
//inline std::basic_istream<CharT,Traits>&
//operator>>(std::basic_istream<CharT,Traits>& in_stream, is_positive_checker<T1> cp)
//{
//  cp.do_check(in_stream);
//  return in_stream;
//}
//
//template <typename T1>
//inline is_positive_checker<T1> assert_positive(T1& var)
//{
//  return is_positive_checker<T1>(var);
//}
//
////--------------------------------------------------------------------------------
//// BINARY PERSISTENCE
////--------------------------------------------------------------------------------
//template <typename It>
//inline void binary_save(std::ostream& out_stream, It begin, It end)
//{
//  typedef typename std::iterator_traits<It>::value_type value_type;
//  size_t size = (size_t) (end - begin);
//  if (size > 0) {
//    char* ptr = (char*) & *begin;
//    out_stream.write(ptr, (std::streamsize) size*sizeof(value_type));
//  }
//}
//
////--------------------------------------------------------------------------------
//template <typename It>
//inline void binary_load(std::istream& in_stream, It begin, It end)
//{
//  typedef typename std::iterator_traits<It>::value_type value_type;
//  size_t size = (size_t) (end - begin);
//  if (size > 0) {
//    char* ptr = (char*) & *begin;
//    in_stream.read(ptr, (std::streamsize) size*sizeof(value_type));
//  }
//}
//
////--------------------------------------------------------------------------------
//template <typename T>
//inline void binary_save(std::ostream& out_stream, const std::vector<T>& v)
//{
//  nta::binary_save(out_stream, v.begin(), v.end());
//}
//
////--------------------------------------------------------------------------------
//template <typename T>
//inline void binary_load(std::istream& in_stream, std::vector<T>& v)
//{
//  nta::binary_load(in_stream, v.begin(), v.end());
//}
//
////--------------------------------------------------------------------------------
//// STL STREAMING OPERATORS
////--------------------------------------------------------------------------------
//
////--------------------------------------------------------------------------------
//// std::pair
////--------------------------------------------------------------------------------
//template <typename T1, typename T2>
//inline std::ostream& operator<<(std::ostream& out_stream, const std::pair<T1, T2>& p)
//{
//  if (io_control.pair_paren)
//    out_stream << "(";
//  out_stream << p.first;
//  out_stream << io_control.pair_sep;
//  out_stream << p.second;
//  if (io_control.pair_paren)
//    out_stream << ")";
//  return out_stream;
//}
//
////--------------------------------------------------------------------------------
//template <typename T1, typename T2>
//inline std::istream& operator>>(std::istream& in_stream, std::pair<T1, T2>& p)
//{
//  in_stream >> p.first >> p.second;
//  return in_stream;
//}
//
////--------------------------------------------------------------------------------
//// std::vector
////--------------------------------------------------------------------------------
//template <typename T, bool>
//struct vector_loader
//{
//  inline void load(size_t, std::istream&, std::vector<T>&);
//};
//
////--------------------------------------------------------------------------------
///**
// * Partial specialization of above functor for primitive types.
// */
//template <typename T>
//struct vector_loader<T, true>
//{
//  inline void load(size_t n, std::istream& in_stream, std::vector<T>& v)
//  {
//    if (io_control.convert_from_sparse == CSR_01) {
//
//      std::fill(v.begin(), v.end(), (T) 0);
//
//      for (size_t i = 0; i != n; ++i) {
//        int index = 0;
//        in_stream >> index;
//        v[index] = (T) 1;
//      }
//
//    } else if (io_control.bit_vector) {
//
//      for (size_t i = 0; i != n; ++i) {
//        float x = 0;
//        in_stream >> x;
//        if (x)
//          v[i] = 1;
//        else
//          v[i] = 0;
//      }
//
//    } else {
//      for (size_t i = 0; i != n; ++i)
//        in_stream >> v[i];
//    }
//  }
//};
//
////--------------------------------------------------------------------------------
///**
// * Partial specialization for non-primitive types.
// */
//template <typename T>
//struct vector_loader<T, false>
//{
//  inline void load(size_t n, std::istream& in_stream, std::vector<T>& v)
//  {
//    for (size_t i = 0; i != n; ++i)
//      in_stream >> v[i];
//  }
//};
//
////--------------------------------------------------------------------------------
///**
// * Factory that will instantiate the right functor to call depending on whether
// * T is a primitive type or not.
// */
//template <typename T>
//inline void vector_load(size_t n, std::istream& in_stream, std::vector<T>& v)
//{
//  vector_loader<T, boost::is_fundamental<T>::value > loader;
//  loader.load(n, in_stream, v);
//}
//
////--------------------------------------------------------------------------------
//template <typename T, bool>
//struct vector_saver
//{
//  inline void save(size_t n, std::ostream& out_stream, const std::vector<T>& v);
//};
//
////--------------------------------------------------------------------------------
///**
// * Partial specialization for primitive types.
// */
//template <typename T>
//struct vector_saver<T, true>
//{
//  inline void save(size_t n, std::ostream& out_stream, const std::vector<T>& v)
//  {
//    if (io_control.output_n_elts)
//      out_stream << n << ' ';
//
//    if (io_control.abbr > 0)
//      n = std::min((size_t) io_control.abbr, n);
//
//    if (io_control.convert_to_sparse == CSR_01) {
//
//      for (size_t i = 0; i != n; ++i)
//        if (!is_zero(v[i]))
//          out_stream << i << ' ';
//
//    } else if (io_control.bit_vector) {
//
//      size_t k = 7;
//      for (size_t i = 0; i != v.size(); ++i) {
//        out_stream << (is_zero(v[i]) ? '0' : '1');
//        if (i == k) {
//          out_stream << ' ';
//          k += 8;
//        }
//      }
//
//    } else {
//
//      for (size_t i = 0; i != n; ++i)
//        out_stream << v[i] << ' ';
//    }
//
//    if (io_control.abbr > 0 && n < v.size()) {
//      size_t rest = v.size() - n;
//      out_stream << "[+" << rest << "/" << count_non_zeros(v) << "]";
//    }
//  }
//};
//
////--------------------------------------------------------------------------------
///**
// * Partial specialization for non-primitive types.
// */
//template <typename T>
//struct vector_saver<T, false>
//{
//  inline void save(size_t n, std::ostream& out_stream, const std::vector<T>& v)
//  {
//    if (io_control.output_n_elts)
//      out_stream << n << ' ';
//
//    if (io_control.abbr > 0)
//      n = std::min((size_t) io_control.abbr, n);
//
//    for (size_t i = 0; i != n; ++i)
//      out_stream << v[i] << ' ';
//
//    if (io_control.abbr > 0 && n < v.size()) {
//      size_t rest = v.size() - n;
//      out_stream << "[+" << rest << "/" << count_non_zeros(v) << "]";
//    }
//  }
//};
//
////--------------------------------------------------------------------------------
///**
// * Factory that will instantiate the right functor to call depending on whether
// * T is a primitive type or not.
// */
//template <typename T>
//inline void vector_save(size_t n, std::ostream& out_stream, const std::vector<T>& v)
//{
//  vector_saver<T, boost::is_fundamental<T>::value> saver;
//  saver.save(n, out_stream, v);
//}
//
////--------------------------------------------------------------------------------
///**
// * Saves the size of the vector.
// */
//template <typename T>
//inline std::ostream& operator<<(std::ostream& out_stream, const std::vector<T>& v)
//{
//  vector_save(v.size(), out_stream, v);
//  return out_stream;
//}
//
////--------------------------------------------------------------------------------
///**
// * Reads in size of the vector, and redimensions it, except if we are reading
// * a sparse binary vector.
// */
//template <typename T>
//inline std::istream&
//operator>>(std::istream& in_stream, std::vector<T>& v)
//{
//  size_t n = 0;
//  in_stream >> n;
//  v.resize(n);
//  vector_load(n, in_stream, v);
//  return in_stream;
//}
//
////--------------------------------------------------------------------------------
///**
// * Doesn't save the size of the buffer itself.
// */
//template <typename T>
//inline std::ostream& operator<<(std::ostream& out_stream, const Buffer<T>& b)
//{
//  vector_save(b.nnz, out_stream, static_cast<const std::vector<T>&>(b));
//  return out_stream;
//}
//
////--------------------------------------------------------------------------------
///**
// * Doesn't set the size of the buffer itself.
// */
//template <typename T>
//inline std::istream& operator>>(std::istream& in_stream, Buffer<T>& b)
//{
//  in_stream >> b.nnz;
//  CPP_UTILS_ASSERT(b.nnz <= b.size());
//  vector_load(b.nnz, in_stream, static_cast<std::vector<T>&>(b));
//  return in_stream;
//}
//
////--------------------------------------------------------------------------------
//// std::set
////--------------------------------------------------------------------------------
//template <typename T1>
//inline std::ostream& operator<<(std::ostream& out_stream, const std::set<T1>& m)
//{
//  typename std::set<T1>::const_iterator
//  it = m.begin(), end = m.end();
//
//  while (it != end) {
//    out_stream << *it << ' ';
//    ++it;
//  }
//
//  return out_stream;
//}
//
////--------------------------------------------------------------------------------
//// std::map
////--------------------------------------------------------------------------------
//template <typename T1, typename T2>
//inline std::ostream& operator<<(std::ostream& out_stream, const std::map<T1, T2>& m)
//{
//  out_stream << m.size() << " ";
//
//  typename std::map<T1, T2>::const_iterator
//  it = m.begin(), end = m.end();
//
//  while (it != end) {
//    out_stream << it->first << ' ' << it->second << ' ';
//    ++it;
//  }
//
//  return out_stream;
//}
//
////--------------------------------------------------------------------------------
//template <typename T1, typename T2>
//inline std::istream& operator>>(std::istream& in_stream, std::map<T1, T2>& m)
//{
//  int size = 0;
//  in_stream >> size;
//
//  for (int i = 0; i != size; ++i) {
//    T1 k; T2 v;
//    in_stream >> k >> v;
//    m.insert(std::make_pair(k, v));
//  }
//
//  return in_stream;
//}
//
////--------------------------------------------------------------------------------
//// MISCELLANEOUS
////--------------------------------------------------------------------------------
//template <typename T>
//inline void show_all_differences(const std::vector<T>& x, const std::vector<T>& y)
//{
//  std::vector<size_t> diffs;
//  find_all_differences(x, y, diffs);
//  std::cout << diffs.size() << " differences: " << std::endl;
//  for (size_t i = 0; i != diffs.size(); ++i)
//    std::cout << "(at:" << diffs[i]
//    << " y=" << x[diffs[i]]
//    << ", ans=" << y[diffs[i]] << ")";
//  std::cout << std::endl;
//}
//
////--------------------------------------------------------------------------------
} // end namespace nta
#endif // CPP_UTILS_STL_IO_HPP
