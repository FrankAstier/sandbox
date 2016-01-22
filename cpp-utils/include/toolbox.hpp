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

#ifndef CPP_UTILS_TOOLBOX_HPP
#define CPP_UTILS_TOOLBOX_HPP

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <iostream>

// Not on win32!!!
#include <sys/time.h>
#include <sys/types.h>
#include <sys/syscall.h>

#include <limits>
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <functional>
#include <numeric>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace utils {

#ifdef WIN32
#undef min
#undef max
#endif

// This include because on darwin86, vDSP provides high quality optimized
// code that exploits SSE.
#ifdef PLATFORM_darwin86
#include <vecLib/vDSP.h>
#endif

//--------------------------------------------------------------------------------
/*
 * Print bits of a char.
 * Today, we can do this and better by creating a std::vector<bool> and sending
 * that to the stream.
 */
void print_binary(std::ostream& stream, const char x) {
  for (int i = 7; i >= 0; --i)
    stream << ((x & (1 << i)) / (1 << i));
}

//--------------------------------------------------------------------------------
/*
 * Print hexadecimal digits of a char.
 * Today, we can do this and better by manipulating the stream with std::hex.
 */
void print_hex(std::ostream& stream, const char x) {
  static const char hexs[] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
                              'A', 'B', 'C', 'D', 'E', 'F'};

  int xx = (int) x;
  for (int i = 16; i >= 1; i /= 16) {
    int p = xx / i;
    xx -= i * p;
    stream << hexs[p];
  }
}

//--------------------------------------------------------------------------------
/**
 * A functions that implements the distance to zero function as a functor.
 * Defining argument_type and result_type directly here instead of inheriting from
 * std::unary_function so that we have an easier time in SWIG Python wrapping.
 */
template<typename T>
struct DistanceToZero {
  typedef T argument_type;
  typedef T result_type;

  inline T operator()(const T& x) const {
    return x >= 0 ? x : -x;
  }
};

//--------------------------------------------------------------------------------
/**
 * A specialization for unsigned ints, where we only need one test (more efficient).
 */
template<>
inline
unsigned int DistanceToZero<unsigned int>::operator()(const unsigned int& x) const { return x; }

//--------------------------------------------------------------------------------
/**
 * Use this functor if T is guaranteed to be positive only.
 */
template<typename T>
struct DistanceToZeroPositive : public std::unary_function<T, T> {
  inline T operator()(const T& x) const {
    return x;
  }
};

//--------------------------------------------------------------------------------
/**
 * This computes the distance to 1 rather than to 0.
 */
template<typename T>
struct DistanceToOne {
  typedef T argument_type;
  typedef T result_type;

  inline T operator()(const T& x) const {
    return x > (T) 1 ? x - (T) 1 : (T) 1 - x;
  }
};

//--------------------------------------------------------------------------------
/**
 * This functor decides whether a number is almost zero or not, using the
 * platform-wide 1e-6.
 */
template<typename D>
struct IsNearlyZero {
  typedef typename D::result_type value_type;

  D dist_;

  // In the case where D::result_type is integral
  // we convert 1e-6 to zero!
  inline IsNearlyZero()
  : dist_() { }

  inline IsNearlyZero(const IsNearlyZero& other)
  : dist_(other.dist_) { }

  inline IsNearlyZero& operator=(const IsNearlyZero& other) {
    if (this != &other)
      dist_ = other.dist_;

    return *this;
  }

  inline bool operator()(const typename D::argument_type& x) const {
    return dist_(x) <= 1e-6;
  }
};

//--------------------------------------------------------------------------------
/**
 * @b Responsibility:
 *  Tell whether an arithmetic value is zero or not, within some precision,
 *  or whether two values are equal or not, within some precision.
 *
 * @b Parameters:
 *  epsilon: accuracy of the comparison
 *
 * @b Returns:
 *  true if |a| <= epsilon
 *  false otherwise
 *
 * @b Requirements:
 *  T arithmetic
 *  T comparable with operator<= AND operator >=
 *
 * @b Restrictions:
 *  Doesn't compile if T is not arithmetic.
 *  In debug mode, asserts if |a| > 10
 *  In debug mode, asserts if a == infinity, quiet_NaN or signaling_NaN
 *
 * @b Notes:
 *  Comparing floating point numbers is a pretty tricky business. Knuth's got
 *  many pages devoted to it in Vol II. Boost::test has a special function to
 *  handle that. One of the problems is that when more bits are allocated to
 *  the integer part as the number gets bigger, there is inherently less
 *  precision in the decimals. But, for comparisons to zero, we can continue
 *  using an absolute epsilon (instead of multiplying epsilon by the number).
 *  In our application, we are anticipating numbers mostly between 0 and 1,
 *  because they represent probabilities.
 */
template<typename T>
inline bool nearlyZero(const T& a, const T& epsilon = T(1e-6)) {
  return a >= -epsilon && a <= epsilon;
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool nearlyEqual(const T& a, const T& b, const T& epsilon = 1e-6) {
  return nearlyZero((b - a), epsilon);
}

//--------------------------------------------------------------------------------
/**
 * A boolean functor that returns true if the element's selected value is found
 * in the (associative) container (needs to support find).
 *
 * Example:
 * =======
 * typedef std::pair<unsigned int, unsigned int> IdxVal;
 * std::list<IdxVal> row;
 * std::set<unsigned int> alreadyGrouped;
 * row.remove_if(SetIncludes<std::set<unsigned int>, select1st<IdxVal> >(alreadyGrouped));
 *
 * Will remove from row (a list of pairs) the pairs whose first element is
 * already contained in alreadyGrouped.
 */
template<typename C1, typename Selector, bool f = false>
struct IsIncluded {
  IsIncluded(const C1& container)
  : sel_(),
    container_(container) { }

  template<typename T>
  inline bool operator()(const T& p) const {
    if (f)
      return container_.find(sel_(p)) == container_.end();
    else
      return container_.find(sel_(p)) != container_.end();
  }

  Selector sel_;
  const C1& container_;
};

//--------------------------------------------------------------------------------
template<class _Pair>
struct select1st : public std::unary_function<_Pair, typename _Pair::first_type> {
  typename _Pair::first_type&
  operator()(_Pair& __x) const { return __x.first; }

  const typename _Pair::first_type&
  operator()(const _Pair& __x) const { return __x.first; }
};

//--------------------------------------------------------------------------------
template<class _Pair>
struct select2nd : public std::unary_function<_Pair, typename _Pair::second_type> {
  typename _Pair::second_type&
  operator()(_Pair& __x) const { return __x.second; }

  const typename _Pair::second_type&
  operator()(const _Pair& __x) const { return __x.second; }
};

//--------------------------------------------------------------------------------
// PAIRS AND TRIPLETS
//--------------------------------------------------------------------------------

/**
 * Lexicographic order:
 * (1,1) < (1,2) < (1,10) < (2,5) < (3,6) < (3,7) ...
 */
template<typename T1, typename T2>
struct lexicographic_2
: public std::binary_function<bool, std::pair<T1, T2>, std::pair<T1, T2> > {
  inline bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
    if (a.first < b.first)
      return true;
    else if (a.first == b.first) if (a.second < b.second)
      return true;
    return false;
  }
};

//--------------------------------------------------------------------------------
/**
 * Order based on the first member of a pair only:
 * (1, 3.5) < (2, 5.6) < (10, 7.1) < (11, 8.5)
 */
template<typename T1, typename T2>
struct less_1st
: public std::binary_function<bool, std::pair<T1, T2>, std::pair<T1, T2> > {
  inline bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
    return a.first < b.first;
  }
};

//--------------------------------------------------------------------------------
/**
 * Order based on the second member of a pair only:
 * (10, 3.5) < (1, 5.6) < (2, 7.1) < (11, 8.5)
 */
template<typename T1, typename T2>
struct less_2nd
: public std::binary_function<bool, std::pair<T1, T2>, std::pair<T1, T2> > {
  inline bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
    return a.second < b.second;
  }
};

//--------------------------------------------------------------------------------
/**
 * Order based on the first member of a pair only:
 * (10, 3.5) > (8, 5.6) > (2, 7.1) > (1, 8.5)
 */
template<typename T1, typename T2>
struct greater_1st
: public std::binary_function<bool, std::pair<T1, T2>, std::pair<T1, T2> > {
  inline bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
    return a.first > b.first;
  }
};

//--------------------------------------------------------------------------------
/**
 * Order based on the second member of a pair only:
 * (10, 3.5) > (1, 5.6) > (2, 7.1) > (11, 8.5)
 */
template<typename T1, typename T2>
struct greater_2nd
: public std::binary_function<bool, std::pair<T1, T2>, std::pair<T1, T2> > {
  inline bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
    return a.second > b.second;
  }
};

//--------------------------------------------------------------------------------
// A class used to work with lists of non-zeros represented in i,j,v format
//--------------------------------------------------------------------------------
/**
 * This class doesn't implement any algorithm, it just stores i,j and v.
 */
template<typename T1, typename T2>
class ijv {
  typedef T1 size_type;
  typedef T2 value_type;

private:
  size_type i_, j_;
  value_type v_;

public:
  inline ijv()
  : i_(0), j_(0), v_(0) { }

  inline ijv(size_type i, size_type j, value_type v)
  : i_(i), j_(j), v_(v) { }

  inline ijv(const ijv& o)
  : i_(o.i_), j_(o.j_), v_(o.v_) { }

  inline ijv& operator=(const ijv& o) {
    i_ = o.i_;
    j_ = o.j_;
    v_ = o.v_;
    return *this;
  }

  inline size_type i() const { return i_; }
  inline size_type j() const { return j_; }
  inline value_type v() const { return v_; }
  inline void i(size_type ii) { i_ = ii; }
  inline void j(size_type jj) { j_ = jj; }
  inline void v(value_type vv) { v_ = vv; }

  //--------------------------------------------------------------------------------
  /**
   * See just above for definition.
   */
  struct lexicographic : public std::binary_function<bool, ijv, ijv> {
    inline bool operator()(const ijv& a, const ijv& b) const {
      if (a.i() < b.i())
        return true;
      else if (a.i() == b.i()) if (a.j() < b.j())
        return true;
      return false;
    }
  };

  //--------------------------------------------------------------------------------
  /**
   * See just above for definition.
   */
  struct less_value : public std::binary_function<bool, ijv, ijv> {
    inline bool operator()(const ijv& a, const ijv& b) const {
      return a.v() < b.v();
    }
  };

  //--------------------------------------------------------------------------------
  /**
   * See just above for definition.
   */
  struct greater_value : public std::binary_function<bool, ijv, ijv> {
    inline bool operator()(const ijv& a, const ijv& b) const {
      return a.v() > b.v();
    }
  };
};

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
//inline size_t unbiased_rand(size_t m)
//{
//  size_t q = m * (numeric_limits<u_int32_t>::max() / m);
//  size_t v;
//  while ((v = arc4random()) > q);
//  return v % m;
//  //return arc4random() % m;
//}
//
////--------------------------------------------------------------------------------
//inline float unbiased_rand_01()
//{
//  return (float)arc4random() / (float)numeric_limits<u_int32_t>::max();
//}

//--------------------------------------------------------------------------------
/**
 * These templates allow the implementation to use the right function depending on
 * the data type, which yield significant speed-ups when working with floats.
 * They also extend STL's arithmetic operations.
 */
//--------------------------------------------------------------------------------
// Unary functions
//--------------------------------------------------------------------------------

template<typename T>
struct Identity : public std::unary_function<T, T> {
  inline T& operator()(T& x) const { return x; }
  inline const T& operator()(const T& x) const { return x; }
};

template<typename T>
struct Negate : public std::unary_function<T, T> {
  inline T operator()(const T& x) const { return -x; }
};

template<typename T>
struct Abs : public std::unary_function<T, T> {
  inline T operator()(const T& x) const { return x > 0.0 ? x : -x; }
};

template<typename T>
struct Square : public std::unary_function<T, T> {
  inline T operator()(const T& x) const { return x * x; }
};

template<typename T>
struct Cube : public std::unary_function<T, T> {
  inline T operator()(const T& x) const { return x * x * x; }
};

template<typename T>
struct Inverse : public std::unary_function<T, T> {
  inline T operator()(const T& x) const { return 1.0 / x; }
};

template<typename T>
struct Sqrt : public std::unary_function<T, T> {
};

template<>
struct Sqrt<float> : public std::unary_function<float, float> {
  inline float operator()(const float& x) const { return sqrtf(x); }
};

template<>
struct Sqrt<double> : public std::unary_function<double, double> {
  inline double operator()(const double& x) const { return sqrt(x); }
};

template<>
struct Sqrt<long double> : public std::unary_function<long double, long double> {
  inline long double operator()(const long double& x) const { return sqrtl(x); }
};

template<typename T>
struct Exp : public std::unary_function<T, T> {
};

template<>
struct Exp<float> : public std::unary_function<float, float> {
  // On x86_64, there is a bug in glibc that makes expf very slow
  // (more than it should be), so we continue using exp on that
  // platform as a workaround.
  // https://bugzilla.redhat.com/show_bug.cgi?id=521190
  // To force the compiler to use exp instead of expf, the return
  // type (and not the argument type!) needs to be double.
  inline float operator()(const float& x) const { return expf(x); }
};

template<>
struct Exp<double> : public std::unary_function<double, double> {
  inline double operator()(const double& x) const { return exp(x); }
};

template<>
struct Exp<long double> : public std::unary_function<long double, long double> {
  inline long double operator()(const long double& x) const { return expl(x); }
};

template<typename T>
struct Log : public std::unary_function<T, T> {
};

template<>
struct Log<float> : public std::unary_function<float, float> {
  inline float operator()(const float& x) const { return logf(x); }
};

template<>
struct Log<double> : public std::unary_function<double, double> {
  inline double operator()(const double& x) const { return log(x); }
};

template<>
struct Log<long double> : public std::unary_function<long double, long double> {
  inline long double operator()(const long double& x) const { return logl(x); }
};

template<typename T>
struct Log2 : public std::unary_function<T, T> {
};

template<>
struct Log2<float> : public std::unary_function<float, float> {
  inline float operator()(const float& x) const {
#ifdef WIN32
    return (float) (log(x) / log(2.0));
#else
    return log2f(x);
#endif
  }
};

template<>
struct Log2<double> : public std::unary_function<double, double> {
  inline double operator()(const double& x) const {
#ifdef WIN32
    return log(x) / log(2.0);
#else
    return log2(x);
#endif
  }
};

template<>
struct Log2<long double> : public std::unary_function<long double, long double> {
  inline long double operator()(const long double& x) const {
#ifdef WIN32
    return log(x) / log(2.0);
#else
    return log2l(x);
#endif
  }
};

template<typename T>
struct Log10 : public std::unary_function<T, T> {
};

template<>
struct Log10<float> : public std::unary_function<float, float> {
  inline float operator()(const float& x) const {
#ifdef WIN32
    return (float) (log(x) / log(10.0));
#else
    return log10f(x);
#endif
  }
};

template<>
struct Log10<double> : public std::unary_function<double, double> {
  inline double operator()(const double& x) const {
#ifdef WIN32
    return log(x) / log(10.0);
#else
    return log10(x);
#endif
  }
};

template<>
struct Log10<long double> : public std::unary_function<long double, long double> {
  inline long double operator()(const long double& x) const {
#ifdef WIN32
    return log(x) / log(10.0);
#else
    return log10l(x);
#endif
  }
};

template<typename T>
struct Log1p : public std::unary_function<T, T> {
};

template<>
struct Log1p<float> : public std::unary_function<float, float> {
  inline float operator()(const float& x) const {
#ifdef WIN32
    return (float) log(1.0 + x);
#else
    return log1pf(x);
#endif
  }
};

template<>
struct Log1p<double> : public std::unary_function<double, double> {
  inline double operator()(const double& x) const {
#ifdef WIN32
    return log(1.0 + x);
#else
    return log1p(x);
#endif
  }
};

template<>
struct Log1p<long double> : public std::unary_function<long double, long double> {
  inline long double operator()(const long double& x) const {
#ifdef WIN32
    return log(1.0 + x);
#else
    return log1pl(x);
#endif
  }
};

/**
 * Numerical approximation of derivative.
 * Error is h^4 y^5/30.
 */
template<typename Float, typename F>
struct Derivative : public std::unary_function<Float, Float> {
  Derivative(const F& f) : f_(f) { }

  F f_;

  /**
   * Approximates the derivative of F at x.
   */
  inline const Float operator()(const Float& x) const {
    const Float h = 1e-6;
    return (-f_(x + 2 * h) + 8 * f_(x + h) - 8 * f_(x - h) + f_(x - 2 * h)) / (12 * h);
  }
};

//--------------------------------------------------------------------------------
// Binary functions
//--------------------------------------------------------------------------------
template<typename T>
struct Assign : public std::binary_function<T, T, T> {
  inline T operator()(T& x, const T& y) const {
    x = y;
    return x;
  }
};

template<typename T>
struct Plus : public std::binary_function<T, T, T> {
  inline T operator()(const T& x, const T& y) const { return x + y; }
};

template<typename T>
struct Minus : public std::binary_function<T, T, T> {
  inline T operator()(const T& x, const T& y) const { return x - y; }
};

template<typename T>
struct Multiplies : public std::binary_function<T, T, T> {
  inline T operator()(const T& x, const T& y) const { return x * y; }
};

template<typename T>
struct Divides : public std::binary_function<T, T, T> {
  inline T operator()(const T& x, const T& y) const { return x / y; }
};

template<typename T>
struct Pow : public std::binary_function<T, T, T> {
};

template<>
struct Pow<float> : public std::binary_function<float, float, float> {
  inline float operator()(const float& x, const float& y) const {
    return powf(x, y);
  }
};

template<>
struct Pow<double> : public std::binary_function<double, double, double> {
  inline double operator()(const double& x, const double& y) const {
    return pow(x, y);
  }
};

template<>
struct Pow<long double>
: public std::binary_function<long double, long double, long double> {
  inline long double operator()(const long double& x, const long double& y) const {
    return powl(x, y);
  }
};

template<typename T>
struct Logk : public std::binary_function<T, T, T> {
};

template<>
struct Logk<float> : public std::binary_function<float, float, float> {
  inline float operator()(const float& x, const float& y) const {
    return logf(x) / logf(y);
  }
};

template<>
struct Logk<double> : public std::binary_function<double, double, double> {
  inline double operator()(const double& x, const double& y) const {
    return log(x) / log(y);
  }
};

template<>
struct Logk<long double>
: public std::binary_function<long double, long double, long double> {
  inline long double operator()(const long double& x, const long double& y) const {
    return logl(x) / logl(y);
  }
};

template<typename T>
struct Max : public std::binary_function<T, T, T> {
  inline T operator()(const T& x, const T& y) const { return x > y ? x : y; }
};

template<typename T>
struct Min : public std::binary_function<T, T, T> {
  inline T operator()(const T& x, const T& y) const { return x < y ? x : y; }
};

/**
 * y = k1*exp(k2*(x-x0)^2) as a functor.
 */
template<typename T>
struct Expk1k2 : public std::unary_function<T, T> {
  T k1_, k2_, x0_;

  inline Expk1k2()
  : k1_(T(1.0 / Sqrt<T>(2 * 3.1415926535))), k2_(-1.0 / .5), x0_(0) { }

  inline Expk1k2(const T& k1, const T& k2, const T& x0)
  : k1_(k1), k2_(k2), x0_(x0) { }

  inline Expk1k2(const Expk1k2& o)
  : k1_(o.k1_), k2_(o.k2_), x0_(o.x0_) { }

  inline Expk1k2& operator=(const Expk1k2& o) {
    if (&o != this) {
      k1_ = o.k1_;
      k2_ = o.k2_;
      x0_ = o.x0_;
    }

    return *this;
  }

  inline T operator()(const T& x) const {
    T v = x - x0_;
    return k1_ * Exp<T>(k2_ * v * v);
  }
};


//--------------------------------------------------------------------------------
/**
 * Gaussian:
 * y = 1/(sigma * sqrt(2*pi)) * exp(-(x-mu)^2/(2*sigma^2)) as a functor.
 */
template<typename T>
struct Gaussian : public std::unary_function<T, T> {
  T k1, k2, mu;

  inline Gaussian(T m, T s)
  : k1(0.0),
    k2(0.0),
    mu(m) {
    // For some reason, SWIG cannot parse 1 / (x), the parentheses in the
    // denominator don't agree with it, so we have to initialize those
    // constants here.
    k1 = 1.0 / sqrt(2.0 * 3.1415926535);
    k2 = -1.0 / (2.0 * s * s);
  }

  inline Gaussian(const Gaussian& o)
  : k1(o.k1), k2(o.k2), mu(o.mu) { }

  inline Gaussian& operator=(const Gaussian& o) {
    if (&o != this) {
      k1 = o.k1;
      k2 = o.k2;
      mu = o.mu;
    }

    return *this;
  }

  inline T operator()(T x) const {
    T v = x - mu;
    return k1 * exp(k2 * v * v);
  }
};

//--------------------------------------------------------------------------------
/**
 * 2D Gaussian
 */
template<typename T>
struct Gaussian2D // : public std::binary_function<T, T, T> (SWIG pb)
{
  T c_x, c_y, s00, s01, s10, s11, s2, k1;

  inline Gaussian2D(T c_x_, T c_y_, T s00_, T s01_, T s10_, T s11_)
  : c_x(c_x_), c_y(c_y_),
    s00(s00_), s01(s01_), s10(s10_), s11(s11_), s2(s10 + s01),
    k1(0.0) {
    // For some reason, SWIG cannot parse 1 / (x), the parentheses in the
    // denominator don't agree with it, so we have to initialize those
    // constants here.
    k1 = 1.0 / (2.0 * 3.1415926535 * sqrt(s00 * s11 - s10 * s01));
    T d = -2.0 * (s00 * s11 - s10 * s01);
    s00 /= d;
    s01 /= d;
    s10 /= d;
    s11 /= d;
    s2 /= d;
  }

  inline Gaussian2D(const Gaussian2D& o)
  : c_x(o.c_x), c_y(o.c_y),
    s00(o.s00), s01(o.s01), s10(o.s10), s11(o.s11), s2(o.s2),
    k1(o.k1) { }

  inline Gaussian2D& operator=(const Gaussian2D& o) {
    if (&o != this) {
      c_x = o.c_x;
      c_y = o.c_y;
      s00 = o.s00;
      s01 = o.s01;
      s10 = o.s10;
      s11 = o.s11;
      s2 = o.s2;
      k1 = o.k1;
    }

    return *this;
  }

  inline T operator()(T x, T y) const {
    T v0 = x - c_x, v1 = y - c_y;
    return k1 * exp(s11 * v0 * v0 + s2 * v0 * v1 + s00 * v1 * v1);
  }
};

////--------------------------------------------------------------------------------
//inline unsigned long long factorial(size_t n) {
//  static boost::unordered_map<size_t, unsigned long long> factorial_cache;
//
//  if (factorial_cache.empty()) {
//
//    unsigned long long fc[] =
//    {0, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800,
//        479001600, 6227020800ull, 87178291200ull};
//
//    for (size_t i = 0; i < 15; ++i)
//      factorial_cache[i] = fc[i];
//  }
//
//  boost::unordered_map<size_t, unsigned long long>::iterator it;
//  it = factorial_cache.find(n);
//  if (it != factorial_cache.end())
//    return it->second;
//
//  unsigned long long f = 1;
//  for (size_t i = 1; i < n; ++i)
//    f *= i;
//  factorial_cache[n] = f;
//  return f;
//}
//
////--------------------------------------------------------------------------------
//inline unsigned long long pow2(size_t n) {
//  static boost::unordered_map<size_t, unsigned long long> pow2_cache;
//
//  if (pow2_cache.empty()) {
//    unsigned long long p2c[] = {1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096,
//                      8192, 16384, 32768, 65536};
//
//    for (size_t i = 0; i <= 16; ++i)
//      pow2_cache[i] = p2c[i];
//  }
//
//  boost::unordered_map<size_t, unsigned long long >::iterator it;
//  it = pow2_cache.find(n);
//  if (it != pow2_cache.end())
//    return it->second;
//
//  unsigned long long val = 65536;
//  for (size_t i = 17; i < n; ++i)
//    val *= 2;
//  pow2_cache[n] = val;
//
//  return val;
//}

//--------------------------------------------------------------------------------
/**
 * Compose two unary functions.
 */
template<typename F1, typename F2>
struct unary_compose
: public std::unary_function<typename F1::argument_type, typename F2::result_type> {
  typedef typename F1::argument_type argument_type;
  typedef typename F2::result_type result_type;

  F1 f1;
  F2 f2;

  inline result_type operator()(const argument_type& x) const {
    return f2(f1(x));
  }
};

//--------------------------------------------------------------------------------
/**
 * Compose an order predicate and a binary selector, so that we can write:
 * sort(x.begin(), x.end(), compose<less<float>, select2nd<pair<int, float> > >());
 * to sort pairs in increasing order of their second element.
 */
template<typename O, typename S>
struct predicate_compose
: public std::binary_function<typename S::argument_type,
typename S::argument_type,
bool> {
  typedef bool result_type;
  typedef typename S::argument_type argument_type;

  O o;
  S s;

  inline result_type operator()(const argument_type& x, const argument_type& y) const {
    return o(s(x), s(y));
  }
};

//--------------------------------------------------------------------------------
/**
 * Having those here allows to hides some STL syntax complexity, and
 * also allows to intercept v for checks (in dividesVal).
 * It also makes writing/reading code easier, not having to deal with
 * the hairy STL type names.
 */
template<typename T>
inline std::binder2nd<Assign<T> > AssignVal(const T& v) {
  return std::bind2nd(Assign<T>(), v);
}

template<typename T>
inline std::binder2nd<Plus<T> > PlusVal(const T& v) {
  return std::bind2nd(Plus<T>(), v);
}

template<typename T>
inline std::binder2nd<Minus<T> > MinusVal(const T& v) {
  return std::bind2nd(Minus<T>(), v);
}

template<typename T>
inline std::binder2nd<Multiplies<T> > MultipliesByVal(const T& v) {
  return std::bind2nd(Multiplies<T>(), v);
}

template<typename T>
inline std::binder2nd<Divides<T> > DividesByVal(const T& v) {
  assert(!nearlyZero(v));

  return std::bind2nd(Divides<T>(), v);
}

template<typename T>
inline std::binder2nd<Pow<T> > PowVal(const T& v) {
  return std::bind2nd(Pow<T>(), v);
}

template<typename T>
inline std::binder2nd<Logk<T> > LogkVal(const T& v) {
  return std::bind2nd(Logk<T>(), v);
}

//--------------------------------------------------------------------------------
/**
 *  When dividing by a value less than min_exponent10, inf will be generated.
 *   numeric_limits<float>::min_exponent10 = -37
 *   numeric_limits<double>::min_exponent10 = -307
 */
template<typename T>
inline bool isSafeForDivision(const T& x) {
  Log <T> log_f;
  return log_f(x) >= std::numeric_limits<T>::min_exponent10;
}

//--------------------------------------------------------------------------------
/**
 * Returns the value passed in or a threshold if the value is >= threshold.
 */
template<typename T>
struct ClipAbove : public std::unary_function<T, T> {
  inline ClipAbove(const T& val)
  : val_(val) { }

  inline ClipAbove(const ClipAbove& c)
  : val_(c.val_) { }

  inline ClipAbove& operator=(const ClipAbove& c) {
    if (this != &c)
      val_ = c.val_;

    return *this;
  }

  inline T operator()(const T& x) const {
    return x >= val_ ? val_ : x;
  }

  T val_;
};

//--------------------------------------------------------------------------------
/**
 * Returns the value passed in or a threshold if the value is < threshold.
 */
template<typename T>
struct ClipBelow : public std::unary_function<T, T> {
  inline ClipBelow(const T& val)
  : val_(val) { }

  inline ClipBelow(const ClipBelow& c)
  : val_(c.val_) { }

  inline ClipBelow& operator=(const ClipBelow& c) {
    if (this != &c)
      val_ = c.val_;

    return *this;
  }

  inline T operator()(const T& x) const {
    return x < val_ ? val_ : x;
  }

  T val_;
};

//--------------------------------------------------------------------------------
typedef unsigned char Byte;

//--------------------------------------------------------------------------------
template<typename T>
inline void print_bits(const T& x) {
  for (int i = sizeof(T) - 1; 0 <= i; --i) {
    unsigned char *b = (unsigned char *) (&x) + i;
    for (int j = 7; 0 <= j; --j)
      std::cout << ((*b & (1 << j)) / (1 << j));
    std::cout << ' ';
  }
}

//--------------------------------------------------------------------------------
// BYTE VECTOR
//--------------------------------------------------------------------------------
/**
 * This is a good compromise between speed and memory for the use cases we have.
 * Going to a real vector of bits is slower when accessing the individual bits,
 * but this vector of bytes can still be fed to the SSE with good results.
 */
struct ByteVector : public std::vector<Byte> {
  inline ByteVector(size_t n = 0)
  : std::vector<Byte>(n, (Byte) 0) { }

  template<typename It>
  inline ByteVector(It begin, It end)
  : std::vector<Byte>((size_t) (end - begin), 0) {
    assert(begin < end);

    for (; begin != end; ++begin)
      (*this)[*begin] = 1;
  }

  /**
   * Use these two functions when converting with a vector of int or float
   * since the byte represenation of the elements in a byte vector is _not_
   * the same as the byte representation of ints and floats.
   */
  template<typename It>
  inline ByteVector(It begin, size_t n)
  : std::vector<Byte>(n, 0) {
    for (size_t i = 0; i != this->size(); ++i)
      (*this)[i] = *begin++ != 0;
  }

  template<typename It>
  inline void toDense(It begin, It end) {
    for (size_t i = 0; i != this->size(); ++i)
      *begin++ = (*this)[i] != 0;
  }
};

//--------------------------------------------------------------------------------
// Direct access with fast erase
//--------------------------------------------------------------------------------
// Records who has been set, so that resetting to zero is fast. Usage pattern
// is to clear the board, do a bunch of sets, look at the board (to test for
// membership for example), then reset the board in the next iteration.
// It trades memory for speed. T is adjustable to be bool (uses vector<bool>,
// 1 bit per element), or int, or ushort, or long, or even float. For
// a membership board (set), on darwin86, unsigned short is fastest on 8/11/2010.
// The clear() method provides a kind of incremental reset.
// Assumes the elements that are set are sparse, that there aren't many compared
// to the size of the board.
template<typename I, typename T>
struct DirectAccess {
  typedef I size_type;
  typedef T value_type;

  std::vector<T> board;
  std::vector<size_type> who;

  inline void resize(size_type m, size_type n = 0) {
    //m = 4 * (m / 4);
    board.resize(m, T());
    who.reserve(n == 0 ? m : n);

    assert(who.size() == 0);
  }

  inline void set(size_type w, const T& v = T(1)) {
    assert(w < board.size());
    assert(v != T());
    assert(who.size() < who.capacity());

    if (board[w] == T()) { // that if doesn't doest much at all (verified)
      who.push_back(w);
      assert(std::set<size_type>(who.begin(), who.end()).size() == who.size());
    }

    board[w] = v;
  }

  inline T get(size_type w) const {
    assert(w < board.size());

    return board[w];
  }

  // Only const operator because the non-const has annoying side-effects
  // that are easily unintended
  inline const T& operator[](size_type w) const {
    return board[w];
  }

  inline void increment(size_type w) {
    assert(w < board.size());

    if (board[w] == T()) {
      who.push_back(w);
      assert(std::set<size_type>(who.begin(), who.end()).size() == who.size());
    }

    ++board[w];
  }

  // If board[w] becomes T() again, we need to update who
  inline void decrement(size_type w) {
    assert(w < board.size());

    if (board[w] == T()) {
      who.push_back(w);
      assert(std::set<size_type>(who.begin(), who.end()).size() == who.size());
    }

    --board[w];

    // To make sure we keep the uniqueness invariant,
    // might be costly if not very sparse?
    if (board[w] == T()) {
      size_type i = 0;
      while (who[i] != w)
        ++i;
      std::swap(who[i], who[who.size() - 1]);
      who.pop_back();
    }
  }

  // v can be anything, < 0, == 0, or > 0
  // If board[w] becomes T() again, we need to update who
  inline void update(size_type w, const value_type& v) {
    assert(w < board.size());

    if (board[w] == T()) {
      who.push_back(w);
      if (std::set<size_type>(who.begin(), who.end()).size() != who.size()) {
//        cout << "Error in update" << endl;
//        cout << board << endl;
//        cout << who << endl;
//        cout << w, v, endl;
        exit(-1);
      }
    }

    board[w] += v;

    // To make sure we keep the uniqueness invariant,
    // might be costly if not very sparse?
    if (board[w] == T()) {
      size_type i = 0;
      while (who[i] != w)
        ++i;
      std::swap(who[i], who[who.size() - 1]);
      who.pop_back();
    }
  }

  // Clear by 4 is a little bit faster, but works only
  // if T() takes exactly 4 bytes.
  inline void clear() {
    size_type *w = &who[0], *w_end = w + who.size();
    //size_type* p = (size_type*) &board[0];
    while (w != w_end)
      //p[*w++] = 0;
      board[*w++] = T();
    who.resize(0);
  }

  // Keep only the value above a certain threshold.
  // Resort the who array optionally.
  // TODO: unit test more
  inline void threshold(const T& t, bool sorted = false) {
    int n = who.size();
    int i = 0;

    while (i < n)
      if (board[who[i]] < t)
        std::swap(who[i], who[--n]);
      else
        ++i;

    who.resize(n);

    if (sorted)
      std::sort(who.begin(), who.end());
  }
};

//--------------------------------------------------------------------------------
// Avoids cost of clearing the board by using multiple colors. Clears only
// every 255 iterations.
// Doesn't keep list of who's on for fast iteration like DirecAccess does.
template<typename I, typename T>
struct Indicator;

template<typename I>
struct Indicator<I, unsigned short> {
  typedef I size_type;

  std::vector<unsigned short> board;
  unsigned short color;

  inline void resize(size_type m) {
    color = 0;
    board.resize(m, color);
  }

  inline void set(size_type w) {
    assert(w < board.size());

    board[w] = color;
  }

  inline bool is_on(size_type w) const {
    assert(w < board.size());

    return board[w] == color;
  }

  inline bool operator[](size_type w) const {
    return is_on(w);
  }

  inline void clear() {
    if (color < std::numeric_limits<unsigned short>::max())
      ++color;
    else {
      color = 0;
      std::fill(board.begin(), board.end(), color);
    }
  }

  template<typename It>
  inline void set_from_sparse(It begin, It end) {
    assert(begin <= end);

    this->clear();
    while (begin != end)
      this->set(*begin++);
  }
};

//--------------------------------------------------------------------------------
/**
 * The first element of each pair is the index of a non-zero, and the second element
 * is the value of the non-zero.
 */
template<typename T1, typename T2>
struct SparseVector : public std::vector<std::pair<T1, T2> > {
  typedef T1 size_type;
  typedef T2 value_type;

  inline SparseVector(size_type s = 0)
  : std::vector<std::pair<T1, T2> >(s) { }
};

//--------------------------------------------------------------------------------
// Unsorted threshold vector
//--------------------------------------------------------------------------------
template<typename T>
inline void unsorted_threshold(std::vector<T>& x, const T& t) {
  int n = x.size();
  int i = 0;

  while (i < n)
    if (x[i] < t)
      std::swap(x[i], x[--n]);
    else
      ++i;

  x.resize(n);
}

//--------------------------------------------------------------------------------
/*
 * Finds simultaneous min and max in less than 2n operations. Best we can do is
 * 3n/2 with the trick below.
 */
template<typename T>
std::pair<size_t, size_t> extrema(const std::vector<T>& x) {
  exit(-1); // buggy, not tested
  if (x.empty())
    return std::make_pair(-1, -1);

  size_t last = 2 * (x.size() / 2);
  size_t argmin = 0, argmax = 0;
  T m = x[0], M = m;

  for (size_t i = 0; i < last; i += 2) {
    if (x[i] < x[i + 1]) { // i+1 runs out of bounds
      if (x[i] < m) {
        m = x[i];
        argmin = i;
      }
      if (M < x[i + 1]) {
        M = x[i + 1];
        argmax = i + 1;
      }
    } else if (x[i + 1] < x[i]) {
      if (M < x[i]) {
        M = x[i];
        argmax = i;
      }
      if (x[i + 1] < m) {
        m = x[i + 1];
        argmin = i + 1;
      }
    }
  }

  if (last < x.size()) {
    if (x[x.size() - 1] < m)
      argmin = x.size() - 1;
    if (M < x[x.size() - 1])
      argmax = x.size() - 1;
  }

  return std::make_pair(argmin, argmax);
}

//--------------------------------------------------------------------------------
// Checks whether the SSE supports the operations we need, i.e. SSE3 and SSE4.
// Returns highest SSE level supported by the CPU: 1, 2, 3 or 41 or 42. It also
// returns -1 if SSE is not present at all.
//
// Refer to Intel manuals for details. Basically, after call to cpuid, the
// interesting bits are set to 1 in either ecx or edx:
// If 25th bit of edx is 1, we have sse: 2^25 = 33554432.
// If 26th bit of edx is 1, we have sse2: 2^26 = 67108864.
// If 0th bit of ecx is 1, we have sse3.
// If 19th bit of ecx is 1, we have sse4.1: 2^19 = 524288.
// If 20th bit of ecx is 1, we have sse4.2: 2^20 = 1048576.
//--------------------------------------------------------------------------------
/*
  static int checkSSE()
  {
  unsigned int f = 1, c = 0, d = 0;

  // #ifdef PLATFORM_win32

  //   __asm {
  //       mov eax, f
  //       cpuid
  //       mov c, ecx
  //       mov d, edx
  //       }

  // #elif defined(PLATFORM_darwin86)

  unsigned int a,b;

  // PIC-compliant asm
  __asm__ __volatile__(
  "pushl %%ebx\n\t"
  "cpuid\n\t"
  "movl %%ebx, %1\n\t"
  "popl %%ebx\n\t"
  : "=a" (a), "=r" (b), "=c" (c), "=d" (d)
  : "a" (f)
  : "cc"
  );
  //#endif

  int ret = -1;
  if (d & 33554432) ret = 1;
  if (d & 67108864) ret = 2;
  if (c & 1) ret = 3;
  if (c & 524288) ret = 41;
  if (c & 1048576) ret = 42;

  return ret;
  }
*/

//--------------------------------------------------------------------------------
// Highest SSE level supported by the CPU: 1, 2, 3 or 41 or 42.
// Note that the asm routines are written for gcc only so far, so we turn them
// off for all platforms except darwin86. Also, they won't work properly on 64 bits
// platforms for now.
//--------------------------------------------------------------------------------
#ifdef PLATFORM_darwin86
static const int SSE_LEVEL = checkSSE();
#else
static const int SSE_LEVEL = -1;
#endif

//--------------------------------------------------------------------------------
// TESTS
//
// TODO: nearly zero for positive numbers
// TODO: is C++ trying to use that for all types??
//--------------------------------------------------------------------------------
template<typename It>
inline bool
nearlyZeroRange(It begin, It end,
                const typename std::iterator_traits<It>::value_type epsilon = 1e-6) {
  {
    assert(begin <= end);
  }

  while (begin != end)
    if (!nearlyZero(*begin++, epsilon))
      return false;
  return true;
}

//--------------------------------------------------------------------------------
template<typename It1, typename It2>
inline bool
nearlyEqualRange(It1 begin1, It1 end1, It2 begin2, It2 end2,
                 const typename std::iterator_traits<It1>::value_type epsilon = 1e-6) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 <= end2 - begin2);
  }

  while (begin1 != end1)
    if (!nearlyEqual(*begin1++, *begin2++, epsilon))
      return false;
  return true;
}

//--------------------------------------------------------------------------------
template<typename Container1, typename Container2>
inline bool
nearlyEqualVector(const Container1& c1, const Container2& c2,
                  const typename Container1::value_type& epsilon = 1e-6) {
  //typedef typename Container1::value_type T1;
  //typedef typename Container2::value_type T2;

  if (c1.size() != c2.size())
    return false;

  return nearlyEqualRange(c1.begin(), c1.end(), c2.begin(), c2.end());
}

//--------------------------------------------------------------------------------
// IS ZERO
//--------------------------------------------------------------------------------
template<typename T>
inline bool is_zero(const T& x) {
  return x == 0;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline bool is_zero(const std::pair<T1, T2>& x) {
  return x.first == 0 && x.second == 0;
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool is_zero(const std::vector<T>& x) {
  for (size_t i = 0; i != x.size(); ++i)
    if (!is_zero(x[i]))
      return false;
  return true;
}

//--------------------------------------------------------------------------------
// DENSE isZero
//--------------------------------------------------------------------------------
/**
 * Scans a binary 0/1 vector to decide whether it is uniformly zero,
 * or if it contains non-zeros (4X faster than C++ loop).
 *
 * If vector x is not aligned on a 16 bytes boundary, the function
 * reverts to slow C++. This can happen when using it with slices of numpy
 * arrays.
 *
 * TODO: find 16 bytes aligned block that can be sent to SSE.
 * TODO: support other platforms than just darwin86 for the fast path.
 * TODO: can we go faster if working on ints rather than floats?
 */
template<typename InputIterator>
inline bool isZero_01(InputIterator x, InputIterator x_end) {
  {
    assert(x <= x_end);
  }

  // On win32, the asm syntax is not correct.
#ifdef PLATFORM_darwin86

  // This test can be moved to compile time using a template with an int
  // parameter, and partial specializations that will match the static
  // const int SSE_LEVEL.
  if (SSE_LEVEL >= 41) { // ptest is a SSE 4.1 instruction

    // n is the total number of floats to process.
    // n1 is the number of floats we can process in parallel using SSE.
    // If x is not aligned on a 4 bytes boundary, we eschew all asm.
    int result = 0;
    int n = (int)(x_end - x);
    int n1 = 0;
    if (((long)x) % 16 == 0)
      n1 = 8 * (n / 8); // we are going to process 2x4 floats at a time

    if (n1 > 0) {

      asm volatile(
                   "pusha\n\t" // save all registers

                   // fill xmm4 with all 1's,
                   // our mask to detect if there are on bits
                   // in the vector or not
                   "subl $16, %%esp\n\t" // allocate 4 floats on the stack
                   "movl $0xffffffff, (%%esp)\n\t" // copy mask 4 times,
                   "movl $0xffffffff, 4(%%esp)\n\t" // then move 16 bytes at once
                   "movl $0xffffffff, 8(%%esp)\n\t" // using movaps
                   "movl $0xffffffff, 12(%%esp)\n\t"
                   "movaps (%%esp), %%xmm4\n\t"
                   "addl $16, %%esp\n\t" // deallocate 4 floats on the stack

                   "0:\n\t"
                   // esi and edi point to the same x, but staggered, so
                   // that we can load 2x4 bytes into xmm0 and xmm1
                   "movaps (%%edi), %%xmm0\n\t" // move 4 floats from x
                   "movaps (%%esi), %%xmm1\n\t" // move another 4 floats from same x
                   "ptest %%xmm4, %%xmm0\n\t"   // ptest first 4 floats, in xmm0
                   "jne 1f\n\t" // jump if ZF = 0, some bit is not zero
                   "ptest %%xmm4, %%xmm1\n\t"   // ptest second 4 floats, in xmm1
                   "jne 1f\n\t" // jump if ZF = 0, some bit is not zero

                   "addl $32, %%edi\n\t"  // jump over 4 floats
                   "addl $32, %%esi\n\t"  // and another 4 floats here
                   "subl $8, %%ecx\n\t" // processed 8 floats
                   "ja 0b\n\t"

                   "movl $0, %0\n\t" // didn't find anything, result = 0 (int)
                   "jmp 2f\n\t" // exit

                   "1:\n\t" // found something
                   "movl $0x1, %0\n\t" // result = 1 (int)

                   "2:\n\t" // exit
                   "popa\n\t" // restore all registers

                   : "=m" (result), "=D" (x)
                   : "D" (x), "S" (x + 4), "c" (n1)
                   :
                   );

      if (result == 1)
        return false;
    }

    // Complete computation by iterating over "stragglers" one by one.
    for (int i = n1; i != n; ++i)
      if (*(x+i) > 0)
        return false;
    return true;

  } else {

    for (; x != x_end; ++x)
      if (*x > 0)
        return false;
    return true;
  }
#else
  for (; x != x_end; ++x)
    if (*x > 0)
      return false;
  return true;
#endif
}

//--------------------------------------------------------------------------------
/**
 * 10X faster than function just above.
 */
inline bool
is_zero_01(const ByteVector& x, size_t begin, size_t end) {
  const Byte *x_beg = &x[begin];
  const Byte *x_end = &x[end];

  // On win32, the asm syntax is not correct.
#ifdef PLATFORM_darwin86

  // This test can be moved to compile time using a template with an int
  // parameter, and partial specializations that will match the static
  // const int SSE_LEVEL.
  if (SSE_LEVEL >= 41) { // ptest is a SSE 4.1 instruction

    // n is the total number of floats to process.
    // n1 is the number of floats we can process in parallel using SSE.
    // If x is not aligned on a 4 bytes boundary, we eschew all asm.
    int result = 0;
    int n = (int)(x_end - x_beg);
    int n1 = 0;
    if (((long)x_beg) % 16 == 0)
      n1 = 32 * (n / 32); // we are going to process 32 bytes at a time

    if (n1 > 0) {

      asm volatile(
                   "pusha\n\t" // save all registers

                   // fill xmm4 with all 1's,
                   // our mask to detect if there are on bits
                   // in the vector or not
                   "subl $16, %%esp\n\t" // allocate 4 floats on the stack
                   "movl $0xffffffff, (%%esp)\n\t" // copy mask 4 times,
                   "movl $0xffffffff, 4(%%esp)\n\t" // then move 16 bytes at once
                   "movl $0xffffffff, 8(%%esp)\n\t" // using movaps
                   "movl $0xffffffff, 12(%%esp)\n\t"
                   "movaps (%%esp), %%xmm4\n\t"
                   "addl $16, %%esp\n\t" // deallocate 4 floats on the stack

                   "0:\n\t"
                   // esi and edi point to the same x, but staggered, so
                   // that we can load 2x4 bytes into xmm0 and xmm1
                   "movaps (%%edi), %%xmm0\n\t" // move 4 floats from x
                   "movaps (%%esi), %%xmm1\n\t" // move another 4 floats from same x
                   "ptest %%xmm4, %%xmm0\n\t"   // ptest first 4 floats, in xmm0
                   "jne 1f\n\t" // jump if ZF = 0, some bit is not zero
                   "ptest %%xmm4, %%xmm1\n\t"   // ptest second 4 floats, in xmm1
                   "jne 1f\n\t" // jump if ZF = 0, some bit is not zero

                   "addl $32, %%edi\n\t"  // jump 32 bytes (16 in xmm0 + 16 in xmm1)
                   "addl $32, %%esi\n\t"  // and another 32 bytes
                   "subl $32, %%ecx\n\t" // processed 32 bytes
                   "ja 0b\n\t"

                   "movl $0, %0\n\t" // didn't find anything, result = 0 (int)
                   "jmp 2f\n\t" // exit

                   "1:\n\t" // found something
                   "movl $0x1, %0\n\t" // result = 1 (int)

                   "2:\n\t" // exit
                   "popa\n\t" // restore all registers

                   : "=m" (result), "=D" (x_beg)
                   : "D" (x_beg), "S" (x_beg + 16), "c" (n1)
                   :
                   );

      if (result == 1)
        return false;
    }

    // Complete computation by iterating over "stragglers" one by one.
    for (int i = n1; i != n; ++i)
      if (*(x_beg+i) > 0)
        return false;
    return true;

  } else {

    for (; x_beg != x_end; ++x_beg)
      if (*x_beg > 0)
        return false;
    return true;
  }
#else
  for (; x_beg != x_end; ++x_beg)
    if (*x_beg > 0)
      return false;
  return true;
#endif
}

//--------------------------------------------------------------------------------
template<typename InIter>
inline bool
positive_less_than(InIter begin, InIter end,
                   const typename std::iterator_traits<InIter>::value_type threshold) {
  {
    assert(begin <= end);
  }

  for (; begin != end; ++begin)
    if (*begin > threshold)
      return false;
  return true;
}

//--------------------------------------------------------------------------------
// N BYTES
//--------------------------------------------------------------------------------
/**
 * For primitive types.
 */
template<typename T>
inline size_t n_bytes(const T&) {
  return sizeof(T);
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline size_t n_bytes(const std::pair<T1, T2>& p) {
  size_t n = n_bytes(p.first) + n_bytes(p.second);
  return n;
}

//--------------------------------------------------------------------------------
/**
 * For more bytes for alignment on x86 with darwin: darwin86 always allocates on
 * 16 bytes boundaries, so the three pointers in the STL vectors (of 32 bits each
 * in -m32), become: 3 * 4 + 4 = 16 bytes. The capacity similarly needs to be
 * adjusted for aligment. On other platforms, the alignment might be different.
 *
 * NOTE/WARNING: this is really "accurate" only on darwin86. And even, it's probably
 * only approximate.
 */
template<typename T>
inline size_t n_bytes(const std::vector<T>& a, size_t alignment = 16) {
  size_t n1 = a.capacity() * sizeof(T);
  if (n1 % alignment != 0)
    n1 = alignment * (n1 / alignment + 1);

  size_t n2 = sizeof(std::vector<T>);
  if (n2 % alignment != 0)
    n2 = alignment * (n2 / alignment + 1);

  return n1 + n2;
}

//--------------------------------------------------------------------------------
template<typename T>
inline float load_factor(const std::vector<T>& x) {
  return (float) x.size() / (float) x.capacity();
}

//--------------------------------------------------------------------------------
template<typename T>
inline void adjust_load_factor(std::vector<T>& x, float target) {
  assert(0.0 <= target && target <= 1.0);

  size_t new_capacity = (size_t) ((float) x.size() / target);

  std::vector<T> y;
  y.reserve(new_capacity);
  y.resize(x.size());
  std::copy(x.begin(), x.end(), y.begin());
  x.swap(y);
}

//--------------------------------------------------------------------------------
// VARIOUS
//--------------------------------------------------------------------------------
//template <typename T>
//inline std::string operator+(const std::string& str, T idx)
//{
//  std::stringstream buff;
//  buff << str << idx;
//  return buff.str();
//}

//--------------------------------------------------------------------------------
template<typename T>
inline void append(const std::vector<T>& a, std::vector<T>& b) {
  b.insert(b.end(), a.begin(), a.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline std::vector<T>& operator+=(std::vector<T>& b, const std::vector<T>& a) {
  append(a, b);
  return b;
}

//--------------------------------------------------------------------------------
template<typename T>
inline void append(const std::set<T>& a, std::set<T>& b) {
  b.insert(a.begin(), a.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline std::set<T>& operator+=(std::set<T>& b, const std::set<T>& a) {
  append(a, b);
  return b;
}

//--------------------------------------------------------------------------------
/*
 * Finds lengths of on and off runs in a vector. Passing ons and offs by reference
 * to avoid copies.
 * If x = 0,1,1,1,0,0,0,1,1,0,1,1,1,1,1,0,0, and threshold == 0,
 * ons = 3,2,5,
 * offs = 1,3,1,2.
 */
template<typename T>
inline void find_runs(const T& threshold,
                      const std::vector<T>& x,
                      std::vector<size_t>& ons,
                      std::vector<size_t>& offs) {
  ons.clear();
  offs.clear();

  if (x.empty())
    return;

  size_t i = 0, j = 0;
  bool on = x[0] > threshold;

  while (i < x.size()) {
    if (on && x[i] <= threshold) {
      on = false;
      ons.push_back(i - j);
      j = i;
    }
    if (!on && x[i] > threshold) {
      on = true;
      offs.push_back(i - j);
      j = i;
    }
    ++i;
  }

  if (i != j) {
    if (on)
      ons.push_back(i - j);
    else
      offs.push_back(i - j);
  }
}

//--------------------------------------------------------------------------------
/*
 * Remove items from list that match the predicate (for example, too close)
 * in when removing, only keep the last one.
 * T is a container type, such as list, that supports size(), erase() and
 * iterators with the same semantics as std::list.
 * If x = 1,2,3, 6,7, 10,11,12, and P returns true if gap >= 2,
 * x becomes: 3, 7, 12.
 */
template<typename T, typename P>
inline void remove_consecutive(T& x, P p) {
  if (x.size() < 2)
    return;

  typename T::iterator it, it2;
  it = x.begin();
  it2 = it;
  ++it2;

  while (it2 != x.end()) {
    if (p(*it2, *it)) { // order of the args to the predicate matters
      it = x.erase(it); // important to grab it
    }
    else
      it = it2;
    ++it2;
  }
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
/*
 * Simple fixed size list, useful for logs of sensor events or user commands,
 * for example. When the maximum allowed size is reached, elements are removed
 * from the front of the log, while new elements are always added at the end
 * (FIFO behavior).
 * TODO: use a vector to save space, and amortize insertions/deletions?
 */
template<typename T>
struct FixedSizeFIFO {
public:
  typedef T value_type;
  typedef size_t size_type;

  typedef typename std::list<T> Container;
  typedef typename Container::iterator iterator;
  typedef typename Container::const_iterator const_iterator;
  typedef typename Container::reverse_iterator reverse_iterator;
  typedef typename Container::const_reverse_iterator const_reverse_iterator;

  FixedSizeFIFO(size_type _maxSize = 10000)
  : mMaxSize(_maxSize),
    mData() { }

  inline bool empty() const { return mData.empty(); }
  inline size_type maxSize() const { return mMaxSize; }
  inline size_type size() const { return mData.size(); }
  inline void clear() { mData.clear(); }

  inline iterator begin() { return mData.begin(); }
  inline iterator end() { return mData.end(); }
  inline const_iterator begin() const { return mData.begin(); }
  inline const_iterator end() const { return mData.end(); }

  inline reverse_iterator rbegin() { return mData.rbegin(); }
  inline reverse_iterator rend() { return mData.rend(); }
  inline const_reverse_iterator rbegin() const { return mData.rbegin(); }
  inline const_reverse_iterator rend() const { return mData.rend(); }

  inline T& back() { return mData.back(); }

  inline void push_back(const T& val) {
    if (mData.size() >= mMaxSize)
      mData.erase(mData.begin());
    mData.push_back(val);
  }

  inline void operator+=(const T& val) {
    push_back(val);
  }

  inline void operator+=(const std::vector<T>& x) {
    for (size_t i = 0; i < x.size(); ++i)
      push_back(x[i]);
  }

  inline iterator erase(iterator it) {
    return mData.erase(it);
  }

  void print(std::ostream& outStream) const {
    outStream << mMaxSize << ' '
    << mData << ' ';
  }

  void read(std::istream& inStream) {
    mData.clear();
    inStream >> mMaxSize
    >> mData;
  }

  size_type clean() {
    size_type n = 0;

    typename std::list<T>::iterator it = mData.begin();
    while (it != mData.end())
      if (!it->invariants()) {
        it = mData.erase(it);
        ++n;
      } else
        ++it;

    return n;
  }

private:
  size_type mMaxSize;
  std::list<T> mData;

  FixedSizeFIFO(const FixedSizeFIFO&);
  FixedSizeFIFO& operator=(const FixedSizeFIFO&);
};

//--------------------------------------------------------------------------------
/*
 * Fixed size log persistence.
 */
template<typename T>
inline std::ostream& operator<<(std::ostream& outStream, const FixedSizeFIFO<T>& fsl) {
  fsl.print(outStream);
  return outStream;
}

//--------------------------------------------------------------------------------
template<typename T>
inline std::istream& operator>>(std::istream& inStream, FixedSizeFIFO<T>& fsl) {
  fsl.read(inStream);
  return inStream;
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
/*
 * A simple kernel smoother where the kernel is Gaussian.
 */
template<typename T>
void smooth(const std::vector<T>& x, std::vector<T>& y, T h) {
  typedef size_t size_type;
  typedef T value_type;

  value_type k2 = -1.0 / (2 * h * h);
  y.resize(x.size(), 0);

  for (size_type i = 0; i < x.size(); ++i) {
    value_type s = 0;
    for (size_type j = 0; j < x.size(); ++j) {
      value_type w = exp(k2 * (j - i) * (j - i));
      y[i] += w * x[j];
      s += w;
    }
    y[i] /= s;
  }
}

//--------------------------------------------------------------------------------
/*
 * Stream only the non-zeros of a dense container.
 */
template<typename T>
inline void save_non_zeros(std::ostream& outStream, const std::vector<T>& x) {
  size_t nnz = count_non_zeros(x);
  outStream << x.size() << ' ' << nnz << ' ';
  for (size_t i = 0; i < x.size(); ++i)
    if (!(x[i] == 0))
      outStream << i << ' ' << x[i] << ' ';
}

//--------------------------------------------------------------------------------
/*
 * Reconstruct a dense container from sparse storage of the non-zeros only.
 */
template<typename T>
inline void read_non_zeros(std::istream& inStream, std::vector<T>& x) {
  x.clear();
  size_t s = 0, nnz = 0;
  inStream >> s >> nnz;
  x.resize(s, 0);
  for (size_t i = 0; i < nnz; ++i) {
    size_t idx;
    T val;
    inStream >> idx >> val;
    x[idx] = val;
  }
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
/*
 * A dense histogram class. Assumes the bins all have the same size.
 * Note that this histogram is not normalized by default.
 */
// class Histogram
// {
// public:
//   typedef size_t size_type;
//   typedef float value_type;

//   Histogram(size_type nBins =0)
//     : bins(nBins, 0)
//   {
//     if (nBins == 0) {
//     }
//   }

//   Histogram(const Histogram& o)
//     : bins(o.bins)
//   {}

//   Histogram& operator=(const Histogram& o)
//   {
//     if (&o != this)
//       bins = o.bins;
//     return *this;
//   }

//   inline size_type size() const { return bins.size(); }

//   inline const std::vector<value_type>& getBins() const
//   {
//     return bins;
//   }

//   /*
//    * Adds the counts contained in binCounts to the bins.
//    */
//   inline void operator+=(const std::vector<value_type>& binCounts)
//   {
//     assert(binCounts.size() == bins.size());

//     for (size_type i = 0; i < bins.size(); ++i)
//       bins[i] += binCounts[i];
//   }

//   /*
//    * Update the given bin with the given value.
//    */
//   inline void update(size_type idx, value_type value)
//   {
//     bins[idx] += value;
//   }

//   /*
//    * Returns the cumulative histogram.
//    */
//   std::vector<value_type> cumulative() const
//   {
//     std::vector<value_type> c(bins.size(), 0);
//     for (size_type i = 1; i < bins.size(); ++i)
//       c[i] += c[i-1];
//     return c;
//   }

//   /*
//    * Persistence.
//    * We save only the non-zeros, because the histograms is
//    * probably sparse, so that saving the non-zeros only is
//    * more economical.
//    */
//   inline void print(std::ostream& outStream) const
//   {
//     save_non_zeros(outStream, bins);
//     outStream << ' ';
//   }

//   inline void read(std::istream& inStream)
//   {
//     read_non_zeros(inStream, bins);
//   }

// private:
//   std::vector<value_type> bins;     // the bins themselves
// };

// //--------------------------------------------------------------------------------
// /*
//  * Histogram persistence.
//  */
// inline std::ostream& operator<<(std::ostream& outStream, const Histogram& hist)
// {
//   hist.print(outStream);
//   return outStream;
// }

// //--------------------------------------------------------------------------------
// inline std::istream& operator>>(std::istream& inStream, Histogram& hist)
// {
//   hist.read(inStream);
//   return inStream;
// }

//--------------------------------------------------------------------------------
// map insert or increment
template<typename T1, typename T2>
inline void increment(std::map<T1, T2>& m, const T1& key, const T2& init = 1) {
  typename std::map<T1, T2>::iterator it = m.find(key);
  if (it != m.end())
    ++it->second;
  else
    m[key] = init;
}

//--------------------------------------------------------------------------------
template<typename K, typename V>
inline bool is_in(const K& key, const std::map<K, V>& m) {
  return m.find(key) != m.end();
}

//--------------------------------------------------------------------------------
// Deriving from std::map to add frequently used functionality
template<typename K, typename V, typename C =std::less<K>,
typename A =std::allocator<std::pair<const K, V> > >
struct dict : public std::map<K, V, C, A> {
  inline bool has_key(const K& key) const {
    return is_in(key, *this);
  }

  // Often useful for histograms, where V is an integral type
  inline void increment(const K& key, const V& init = 1) {
    increment(*this, key, init);
  }

  // Inserts once in the map, or return false if already inserted
  // (saves having to write find(...) == this->end())
  inline bool insert_once(const K& key, const V& v) {
    if (has_key(key))
      return false;
    else
      this->insert(std::make_pair(key, v));
    return true;
  }

  /*
  // Returns an existing value for the key, if it is in the dict already,
  // or creates one and returns it. (operator[] on std::map does that?)
  inline V& operator(const K& key)
  {
  iterator it = this->find(key);
  if (key == end) {
  (*this)[key] = V();
  return (*this)[key];
  } else
  return *it;
  }
  */
};

//--------------------------------------------------------------------------------
// INIT LIST - it's in boost assignment
//--------------------------------------------------------------------------------
// template <typename T>
// struct init_list
// {
//   T& v;

//   inline init_list(T& v_ref) : v(v_ref) {}
//   inline init_list(const init_list& o) : v(o.v) {}

//   inline init_list& operator=(const init_list& o)
//   { v(o.v); return *this; }

//   template <typename T2>
//   inline init_list<T>& operator,(const T2& x)
//   {
//     v.push_back(x);
//     return *this;
//   }
// };

// //--------------------------------------------------------------------------------
// template <typename T, typename T2>
// inline init_list<T> operator+=(T& v, const T2& x)
// {
//   v.push_back(x);
//   return init_list<T>(v);
// }

// //--------------------------------------------------------------------------------
// // TODO: merge with preceding by changing parametrization?
// //--------------------------------------------------------------------------------
// template <typename T>
// struct set_init_list
// {
//   std::set<T>& v;

//   inline set_init_list(std::set<T>& v_ref) : v(v_ref) {}
//   inline set_init_list(const set_init_list& o) : v(o.v) {}

//   inline set_init_list& operator=(const set_init_list& o)
//   { v(o.v); return *this; }

//   template <typename T2>
//   inline set_init_list<T>& operator,(const T2& x)
//   {
//     v.insert((T)x);
//     return *this;
//   }
// };

// //--------------------------------------------------------------------------------
// template <typename T, typename T2>
// inline set_init_list<T> operator+=(std::set<T>& v, const T2& x)
// {
//   v.insert((T)x);
//   return set_init_list<T>(v);
// }

// //--------------------------------------------------------------------------------
// /*
//  * A class to initialize lists conveniently.
//  */
// template <typename T>
// struct list_init_list
// {
//   std::list<T>& v;

//   inline list_init_list(std::list<T>& v_ref) : v(v_ref) {}
//   inline list_init_list(const list_init_list& o) : v(o.v) {}

//   inline list_init_list& operator=(const list_init_list& o)
//   { v(o.v); return *this; }

//   template <typename T2>
//   inline list_init_list<T>& operator,(const T2& x)
//   {
//     v.push_back((T)x);
//     return *this;
//   }
// };

// //--------------------------------------------------------------------------------
// /*
//  * Works with list_init_list to initialize a list conveniently.
//  */
// template <typename T, typename T2>
// inline list_init_list<T> operator+=(std::list<T>& v, const T2& x)
// {
//   v.push_back((T)x);
//   return list_init_list<T>(v);
// }

//--------------------------------------------------------------------------------
// FIND IN VECTOR
//--------------------------------------------------------------------------------
// T1 and T2 to get around constness with pointers
template<typename T1, typename T2>
inline int find_index(const T1& x, const std::vector<T2>& v) {
  for (size_t i = 0; i != v.size(); ++i)
    if (v[i] == x)
      return (int) i;
  return -1;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline int find_index(const T1& x, const std::vector<std::pair<T1, T2> >& v) {
  for (size_t i = 0; i != v.size(); ++i)
    if (v[i].first == x)
      return (int) i;
  return -1;
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool not_in(const T& x, const std::vector<T>& v) {
  return std::find(v.begin(), v.end(), x) == v.end();
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline bool not_in(const T1& x, const std::vector<std::pair<T1, T2> >& v) {
  typename std::vector<std::pair<T1, T2> >::const_iterator it;
  for (it = v.begin(); it != v.end(); ++it)
    if (it->first == x)
      return false;
  return true;
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool not_in(const T& x, const std::set<T>& s) {
  return s.find(x) == s.end();
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool is_in(const T& x, const std::vector<T>& v) {
  return !not_in(x, v);
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline bool is_in(const T1& x, const std::vector<std::pair<T1, T2> >& v) {
  return !not_in(x, v);
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool is_in(const T& x, const std::set<T>& s) {
  return !not_in(x, s);
}

//--------------------------------------------------------------------------------
/*
 * Checks whether a range is sorted or not.
 */
template<typename It>
inline bool is_sorted(It begin, It end, bool ascending = true, bool unique = true) {
  for (It prev = begin, it = ++begin; it < end; ++it, ++prev)

    if (ascending) {
      if (unique) {
        if (!(*prev < *it))
          return false;
      } else {
        if (*it < *prev)
          return false;
      }
    } else {
      if (unique) {
        if (!(*it < *prev))
          return false;
      } else {
        if (*prev < *it)
          return false;
      }
    }

  return true;
}

//--------------------------------------------------------------------------------
/*
 * Checks whether a whole container is sorted or not.
 */
template<typename C>
inline bool is_sorted(const C& x, bool ascending = true, bool unique = true) {
  return is_sorted(x.begin(), x.end(), ascending, unique);
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool operator==(const std::vector<T>& a, const std::vector<T>& b) {
  if (a.size() != b.size())
    return false;
  for (size_t i = 0; i != a.size(); ++i)
    if (a[i] != b[i])
      return false;
  return true;
}

//--------------------------------------------------------------------------------
template<typename T>
inline bool operator!=(const std::vector<T>& a, const std::vector<T>& b) {
  return !(a == b);
}

//--------------------------------------------------------------------------------
/*
 * A runlength-encoded vector.
 */
template<typename T>
struct RLEVector : public std::vector<std::pair<T, size_t> > {
};

//--------------------------------------------------------------------------------
/*
 * Compresses a vector. For example:
 * x = 1 1 1 1 2 2 2 3 3 3 3 3 <=> y = (1,4) (2,3) (3,5)
 */
template<typename T>
inline void
RLECompressVector(const std::vector<T>& x, RLEVector<T>& y) {
  y.clear();
  size_t i = 0;
  while (i < x.size()) {
    T val = x[i];
    size_t counter = 0;
    while (i < x.size() && x[i] == val) {
      ++counter;
      ++i;
    }
    y.push_back(std::make_pair(val, counter));
  }
}

//--------------------------------------------------------------------------------
/*
 * Inverse of RLECompressVector.
 */
template<typename T>
inline void
RLEUncompressVector(const RLEVector<T>& y, std::vector<T>& x) {
  x.clear();
  for (size_t i = 0; i < y.size(); ++i) {
    T val = y[i].first;
    size_t counter = y[i].second;
    for (size_t j = 0; j < counter; ++j)
      x.push_back(val);
  }
}

//--------------------------------------------------------------------------------
/**
 * Proxy for an insert iterator that allows inserting at the second element
 * when iterating over a container of pairs.
 */
template<typename Iterator>
struct inserter_second {
  typedef typename std::iterator_traits<Iterator>::value_type pair_type;
  typedef typename pair_type::second_type second_type;
  typedef second_type value_type;

  Iterator it;

  inline inserter_second(Iterator _it) : it(_it) { }
  inline second_type& operator*() { return it->second; }
  inline void operator++() { ++it; }
};

template<typename Iterator>
inserter_second<Iterator> insert_2nd(Iterator it) {
  return inserter_second<Iterator>(it);
}

//--------------------------------------------------------------------------------
/**
 * Proxy for an insert iterator that allows inserting at the second element when
 * iterating over a container of pairs, while setting the first element to the
 * current index value (watch out if iterator passed to constructor is not
 * pointing to the beginning of the container!)
 */
template<typename Iterator>
struct inserter_second_incrementer_first {
  typedef typename std::iterator_traits<Iterator>::value_type pair_type;
  typedef typename pair_type::second_type second_type;
  typedef second_type value_type;

  Iterator it;
  size_t i;

  inline inserter_second_incrementer_first(Iterator _it)
  : it(_it), i(0) { }
  inline second_type& operator*() { return it->second; }
  inline void operator++() {
    it->first = i++;
    ++it;
  }
};

template<typename Iterator>
inserter_second_incrementer_first<Iterator> insert_2nd_inc(Iterator it) {
  return inserter_second_incrementer_first<Iterator>(it);
}

//--------------------------------------------------------------------------------
/**
 * Sparse dot product: x and y are vectors that contain the indices of non-zeros,
 * and this functions counts the number of matches between those two vectors.
 */
template<typename It1, typename T2>
inline size_t sparse_dot(It1 x_begin, It1 x_end, const std::vector<T2>& y) {
  size_t n2 = y.size(), i2 = 0;
  size_t s = 0;

  while (x_begin != x_end && i2 != n2)
    if (*x_begin < y[i2]) {
      ++x_begin;
    } else if (y[i2] < *x_begin) {
      ++i2;
    } else {
      ++s;
      ++x_begin;
      ++i2;
    }

  return s;
}

//--------------------------------------------------------------------------------
/**
 * Sparse dot product: x and y are vectors that contain the indices of non-zeros,
 * and this functions counts the number of matches between those two vectors.
 */
template<typename T1, typename T2>
inline T2 sparse_dot(const std::vector<T1>& x, const std::vector<T2>& y) {
  size_t n1 = x.size(), n2 = y.size(), i1 = 0, i2 = 0;
  T2 s = 0;

  while (i1 != n1 && i2 != n2)
    if (x[i1] < y[i2]) {
      ++i1;
    } else if (y[i2] < x[i1]) {
      ++i2;
    } else {
      ++s;
      ++i1;
      ++i2;
    }

  return s;
}

//--------------------------------------------------------------------------------
inline float dot(const float *x, const float *x_end, const float *y) {
  float result = 0;
#ifdef PLATFORM_darwin86
  vDSP_dotpr(x, 1, y, 1, &result, (x_end - x));
#else
  for (; x != x_end; ++x, ++y)
    result += *x * *y;
#endif
  return result;
}

//--------------------------------------------------------------------------------
// copy
//--------------------------------------------------------------------------------
template<typename It1, typename It2>
inline void copy(It1 begin, It1 end, It2 out_begin, It2 out_end) {
  std::copy(begin, end, out_begin);
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline void copy(const T1& a, T2& b) {
  b.resize(a.size());
  copy(a.begin(), a.end(), b.begin(), b.end());
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline void copy(const std::vector<T1>& a, size_t n, std::vector<T2>& b, size_t o) {
  assert(o + n <= b.size());
  std::copy(a.begin(), a.begin() + n, b.begin() + o);
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline void copy(const std::vector<T1>& a, size_t i, size_t j, std::vector<T2>& b) {
  std::copy(a.begin() + i, a.begin() + j, b.begin() + i);
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline void copy(const std::vector<T1>& a, std::vector<T2>& b, size_t offset) {
  assert(offset + a.size() <= b.size());
  std::copy(a.begin(), a.end(), b.begin() + offset);
}

//--------------------------------------------------------------------------------
template<typename I, typename T>
inline void copy_indices(const SparseVector<I, T>& x, std::vector<I>& y) {
  assert(x.size() <= y.size());

  for (size_t i = 0; i != x.size(); ++i)
    y[i] = x[i].first;
  y.size() = x.size();
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2, typename T3>
inline void copy_first(const std::vector<std::pair<T1, T2> >& x, std::vector<T3>& y) {
  y.resize(x.size());
  for (size_t i = 0; i != x.size(); ++i)
    y[i] = x[i].first;
}

//--------------------------------------------------------------------------------
// TO DENSE
//--------------------------------------------------------------------------------
template<typename It1, typename It2>
inline void to_dense_01(It1 ind, It1 ind_end, It2 dense, It2 dense_end) {
  {
    assert(ind <= ind_end);
    assert(dense <= dense_end);
    assert(ind_end - ind <= dense_end - dense);
  }

  typedef typename std::iterator_traits<It2>::value_type value_type;

  // TODO: make faster with single pass?
  // (but if's for all the elements might be slower)
  std::fill(dense, dense_end, (value_type) 0);

  for (; ind != ind_end; ++ind)
    *(dense + *ind) = (value_type) 1;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline void to_dense_01(const std::vector<T1>& sparse, std::vector<T2>& dense) {
  to_dense(sparse.begin(), sparse.end(), dense.begin(), dense.end());
}

//--------------------------------------------------------------------------------
template<typename It, typename T>
inline void to_dense_01(It begin, It end, std::vector<T>& dense) {
  to_dense_01(begin, end, dense.begin(), dense.end());
}

//--------------------------------------------------------------------------------
template<typename T, typename OutIt>
inline void to_dense_01(const std::vector<T>& buffer, OutIt y, OutIt y_end) {
  typedef typename std::iterator_traits<OutIt>::value_type value_type;

  std::fill(y, y_end, (value_type) 0);

  for (size_t i = 0; i != buffer.size(); ++i)
    y[buffer[i]] = (value_type) 1;
}

//--------------------------------------------------------------------------------
template<typename I, typename T, typename OutIt>
inline void to_dense_1st_01(const SparseVector<I, T>& x, OutIt y, OutIt y_end) {
  typedef typename std::iterator_traits<OutIt>::value_type value_type;

  std::fill(y, y_end, (value_type) 0);

  for (size_t i = 0; i != x.size(); ++i)
    y[x[i].first] = (value_type) 1;
}

//--------------------------------------------------------------------------------
/**
 * Converts a sparse range described with indices and values to a dense
 * range.
 */
template<typename It1, typename It2, typename It3>
inline void to_dense(It1 ind, It1 ind_end, It2 nz, It2 nz_end,
                     It3 dense, It3 dense_end) {
  {
    assert(ind <= ind_end);
    assert(dense <= dense_end);
    assert(ind_end - ind <= dense_end - dense);
    assert(nz_end - nz == ind_end - ind);
  }

  typedef typename std::iterator_traits<It3>::value_type value_type;

  std::fill(dense, dense + (ind_end - ind), (value_type) 0);

  for (; ind != ind_end; ++ind, ++nz)
    *(dense + *ind) = *nz;
}

//--------------------------------------------------------------------------------
/**
 * Needs non-zeros to be sorted!
 */
template<typename It>
inline void in_place_sparse_to_dense_01(int n, It begin, It end) {
  for (int i = n - 1; i >= 0; --i) {
    int p = (int) *(begin + i);
    std::fill(begin + p, end, 0);
    *(begin + p) = 1;
    end = begin + p;
  }

  std::fill(begin, end, 0);
}

//--------------------------------------------------------------------------------
/**
 * Pb with size of the vectors?
 * n is the number of non-zeros stored originally in x.
 */
template<typename T>
inline void in_place_sparse_to_dense_01(int n, std::vector<T>& x) {
  in_place_sparse_to_dense_01(n, x.begin(), x.end());
}

//--------------------------------------------------------------------------------
/**
 * Converts a sparse range stored in a dense vector into an (index,value)
 * representation.
 *
 * @param begin
 * @param end
 * @param ind
 * @param nz
 * @param eps
 */
template<typename It1, typename It2, typename It3>
inline void
from_dense(It1 begin, It1 end, It2 ind, It3 nz,
           typename std::iterator_traits<It1>::value_type eps = 1e-6) {
  {
    assert(begin <= end);
  }

  typedef size_t size_type;
  typedef typename std::iterator_traits<It1>::value_type value_type;

  Abs <value_type> abs_f;

  for (It1 it = begin; it != end; ++it) {
    value_type val = *it;
    if (abs_f(val) > eps) {
      *ind = (size_type) (it - begin);
      *nz = val;
      ++ind;
      ++nz;
    }
  }
}

//--------------------------------------------------------------------------------
/**
 * Converts a sparse vector represented in extenso to a sparse vector that
 * contains only the indices of the non-zeros (ignores the values of the non-zeros).
 */
template<typename T>
inline void in_place_dense_to_sparse_01(std::vector<T>& x, bool adjust = true) {
  typename std::vector<T>::iterator begin = x.begin(), end = x.end();

  size_t last_sparse = 0;
  for (size_t i = 0; i != x.size(); ++i)
    if (x[i] != 0)
      x[last_sparse++] = i;

  if (adjust)
    x.resize(last_sparse);
  else
    std::fill(x.begin() + last_sparse, x.end(), 0);
}

//--------------------------------------------------------------------------------
// erase from vector
//--------------------------------------------------------------------------------
template<typename T>
inline void moveToBackAndPop(size_t idx, std::vector<T>& x) {
  std::swap(idx, x[x.size() - 1]);
  x.pop_back();
}

//--------------------------------------------------------------------------------
/**
 * Erases a value from a vector.
 *
 * The STL process to really remove a value from a vector is tricky.
 *
 * @param v the vector
 * @param val the value to remove
 */
template<typename T>
inline void remove(const T& del, std::vector<T>& v) {
  v.erase(std::remove(v.begin(), v.end(), del), v.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline void remove(const std::vector<T>& del, std::vector<T>& b) {
  for (size_t i = 0; i != del.size(); ++i)
    remove(del[i], b);
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline void remove_pair(const T1& key, std::vector<std::pair<T1, T2> >& v) {
  typename std::vector<std::pair<T1, T2> >::const_iterator it;
  for (it = v.begin(); it != v.end() && it->first != key; ++it);
  remove(*it, v);
}

//--------------------------------------------------------------------------------
template<typename T>
inline void remove_from_end(const T& elt, std::vector<T>& a) {
  for (int i = a.size() - 1; i >= 0; --i) {
    if (a[i] == elt) {
      for (size_t j = i; j < a.size() - 1; ++j)
        a[j] = a[j + 1];
      a.resize(a.size() - 1);
      return;
    }
  }
}

//--------------------------------------------------------------------------------
/**
 * Remove several elements from a vector, the elements to remove being specified
 * by their index (in del). After this call, a's size is reduced. Requires
 * default constructor on T to be defined (for resize). O(n).
 */
template<typename I, typename T>
inline void
remove_at(const std::vector<I>& del, std::vector<T>& a) {
  assert(std::set<I>(del.begin(), del.end()).size() == del.size());

  if (del.empty())
    return;

  size_t old = del[0] + 1, cur = del[0], d = 1;

  while (old < a.size() && d < del.size()) {
    if (old == (size_t) del[d]) {
      ++d;
      ++old;
    } else if ((size_t) del[d] < old) {
      ++d;
    } else {
      a[cur++] = a[old++];
    }
  }

  while (old < a.size())
    a[cur++] = a[old++];

  a.resize(a.size() - del.size());
}

//--------------------------------------------------------------------------------
/**
 * Given a vector of indices, removes the elements of a at those indices (indices
 * before any removal is carried out), where a is a vector of pairs.
 *
 * Need to pass in non-empty vector of sorted, unique indices to delete.
 */
template<typename I, typename T1, typename T2>
inline void
remove_at(const std::vector<I>& del, std::vector<std::pair<T1, T2> >& a) {
  assert(std::set<I>(del.begin(), del.end()).size() == del.size());

  if (del.empty())
    return;

  size_t old = del[0] + 1, cur = del[0], d = 1;

  while (old < a.size() && d < del.size()) {
    if (old == (size_t) del[d]) {
      ++d;
      ++old;
    } else if ((size_t) del[d] < old) {
      ++d;
    } else {
      a[cur++] = a[old++];
    }
  }

  while (old < a.size())
    a[cur++] = a[old++];

  a.resize(a.size() - del.size());
}

//--------------------------------------------------------------------------------
/**
 * Finds index of elt in ref, and removes corresponding element of a.
 */
template<typename T1, typename T2>
inline void remove(const T2& elt, std::vector<T1>& a, const std::vector<T2>& ref) {
  a.erase(a.begin() + find_index(elt, ref));
}
//--------------------------------------------------------------------------------
template<typename T>
inline void remove(const std::vector<T>& del, std::set<T>& a) {
  for (size_t i = 0; i != del.size(); ++i)
    a.erase(del[i]);
}

//--------------------------------------------------------------------------------
template<typename T>
inline void remove(const std::set<T>& y, std::vector<T>& x) {
  std::vector<T> del;

  for (size_t i = 0; i != x.size(); ++i)
    if (y.find(x[i]) != y.end()) {
      assert(not_in(x[i], del));
      del.push_back(x[i]);
    }

  remove(del, x);
}

//--------------------------------------------------------------------------------
template<typename T>
inline std::set<T>& operator-=(std::set<T>& a, const std::vector<T>& b) {
  remove(b, a);
  return a;
}

//--------------------------------------------------------------------------------
template<typename T>
inline std::vector<T>& operator-=(std::vector<T>& a, const std::vector<T>& b) {
  remove(b, a);
  return a;
}

//--------------------------------------------------------------------------------
template<typename T>
inline std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b) {
  std::vector<T> r(a.size(), 0);
  for (size_t i = 0; i != a.size(); ++i)
    r[i] = a[i] - b[i];
  return r;
}

//--------------------------------------------------------------------------------
// DIFFERENCES
//--------------------------------------------------------------------------------
/**
 * Returns a vector that contains the indices of the positions where x and y
 * have different values.
 */
template<typename T>
inline void
find_all_differences(const std::vector<T>& x, const std::vector<T>& y,
                     std::vector<size_t>& diffs) {
  assert(x.size() == y.size());
  diffs.clear();
  for (size_t i = 0; i != x.size(); ++i)
    if (x[i] != y[i])
      diffs.push_back(i);
}

//--------------------------------------------------------------------------------
// fill
//--------------------------------------------------------------------------------
/**
 * Fills a container with the given value.
 *
 * @param a
 * @param val
 */
template<typename T>
inline void fill(T& a, const typename T::value_type& val) {
  std::fill(a.begin(), a.end(), val);
}

//--------------------------------------------------------------------------------
/**
 * Zeroes out a range.
 *
 * @param begin
 * @param end
 */
template<typename It>
inline void zero(It begin, It end) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type T;

  for (; begin != end; ++begin)
    *begin = T(0);
}

//--------------------------------------------------------------------------------
/**
 * Zeroes out a whole container.
 *
 * @param a the container
 */
template<typename T>
inline void zero(T& a) {
  zero(a.begin(), a.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline void set_to_zero(T& a) {
  zero(a);
}

//--------------------------------------------------------------------------------
template<typename T>
inline void set_to_zero(std::vector<T>& a, size_t begin, size_t end) {
  zero(a.begin() + begin, a.begin() + end);
}

//--------------------------------------------------------------------------------
/**
 * Fills a range with ones.
 *
 * @param begin
 * @param end
 */
template<typename It>
inline void ones(It begin, It end) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type T;

  for (; begin != end; ++begin)
    *begin = T(1);
}

//--------------------------------------------------------------------------------
/**
 * Fills a container with ones.
 *
 * @param a the container
 */
template<typename T>
inline void ones(T& a) {
  ones(a.begin(), a.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline void set_to_one(std::vector<T>& a) {
  ones(a);
}

//--------------------------------------------------------------------------------
template<typename T>
inline void set_to_one(std::vector<T>& a, size_t begin, size_t end) {
  ones(a.begin() + begin, a.begin() + end);
}

//--------------------------------------------------------------------------------
/**
 * Fills the container with a range of values.
 */
template<typename T>
inline void generate_range(T& t,
                           typename T::value_type start,
                           typename T::value_type end,
                           typename T::value_type increment = 1) {
  std::insert_iterator<T> it(t, t.begin());

  for (typename T::value_type i = start; i < end; i += increment, ++it)
    *it = i;
}

//--------------------------------------------------------------------------------
inline std::vector<size_t> int_range(size_t n) {
  std::vector<size_t> p(n);
  for (size_t i = 0; i != n; ++i)
    p[i] = i;
  return p;
}

//--------------------------------------------------------------------------------
/**
 * Initializes a range with the uniform distribution.
 *
 * @param begin beginning of the range
 * @param end one past the end of the range
 * @param val the value to which the sum of the range will be equal to
 */
template<typename It>
inline void
uniform_range(It begin, It end,
              typename std::iterator_traits<It>::value_type val = 1) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  std::fill(begin, end, (value_type) 1);
  normalize(begin, end, val);
}

//--------------------------------------------------------------------------------
/**
 * Initializes a container with the uniform distribution.
 *
 * @param a the container
 * @param val the value for normalization
 */
template<typename C>
inline void uniform_range(C& a, typename C::value_type val = 1) {
  uniform_range(a.begin(), a.end(), val);
}

//--------------------------------------------------------------------------------
/**
 * Sets a range to 0, except for a single value at pos, which will be equal to val.
 *
 * @param pos the position of the single non-zero value
 * @param begin
 * @param end
 * @param val the value of the non-zero value in the range
 */
template<typename It>
inline void dirac(size_t pos, It begin, It end,
                  typename std::iterator_traits<It>::value_type val = 1) {
  {
    assert(begin <= end);
    assert(0 <= pos && pos < (size_t) (end - begin));
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  std::fill(begin, end, (value_type) 0);
  *(begin + pos) = val;
}

//--------------------------------------------------------------------------------
/**
 * Sets a range to 0, except for a single value at pos, which will be equal to val.
 *
 * @param pos the position of the single non-zero value
 * @param c the container
 * @param val the value of the Dirac
 */
template<typename C>
inline void dirac(size_t pos, C& c, typename C::value_type val = 1) {
  {
    assert(pos >= 0 && pos < c.size());
  }

  dirac(pos, c.begin(), c.end(), val);
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// STATISTICS
//--------------------------------------------------------------------------------
// Computes cdf for arbitrary pdf passed as functor
//--------------------------------------------------------------------------------
template<typename T, typename F>
inline void cumulative(T a, T b, std::vector<T>& cdf, F pdf) {
  {
    assert(a < b);
    assert(!cdf.empty());
    assert(0 < double(b - a) / cdf.size());
  }

  size_t n = cdf.size(), i = 0;
  double step = (b - a) / n;
  for (double x = a; x <= b; x += step, ++i) {
    cdf[i] = pdf(x);
    assert(0 <= cdf[i] && cdf[i] <= 1);
  }
  for (size_t i = 1; i != n; ++i)
    cdf[i] += cdf[i - 1];
  for (size_t i = 0; i != n; ++i)
    cdf[i] /= cdf[n - 1];
  cdf[n - 1] = 1.0;

  assert(0 <= cdf[0]);
  for (size_t i = 1; i != n; ++i) {
    assert(0 <= cdf[i]);
    assert(cdf[i - 1] <= cdf[i]);
  }
}

//--------------------------------------------------------------------------------
/**
 * Computes the CDF of the given range seen as a discrete PMF.
 *
 * @param begin1 the beginning of the discrete PMF range
 * @param end1 one past the end of the discrete PMF range
 * @param begin2 the beginning of the CDF range
 */
template<typename It1, typename It2>
inline void cumulative(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typedef typename std::iterator_traits<It2>::value_type value_type;

  It2 prev = begin2;
  *begin2++ = (value_type) *begin1++;
  for (; begin1 < end1; ++begin1, ++begin2, ++prev)
    *begin2 = *prev + (value_type) *begin1;
}

//--------------------------------------------------------------------------------
/**
 * Computes the CDF of a discrete PMF.
 *
 * @param pmf the PMF
 * @param cdf the CDF
 */
template<typename C1, typename C2>
inline void cumulative(const C1& pmf, C2& cdf) {
  cumulative(pmf.begin(), pmf.end(), cdf.begin(), cdf.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline void compute_erf(T mu, T sigma, T a, T b, std::vector<T>& cdf) {
  cumulative(a, b, cdf, Gaussian<T>(mu, sigma));
}

//--------------------------------------------------------------------------------
/*
 * From a monotonic function discretized in array f, computes the inverse into
 * array inv_f. The inverse can be discretized if T2 is integral.
 * Uses array f to interpolate the original function when needed.
 * The method is not very precise, nor very fast!
 * Useful to sample from distributions, where the cdf can be computed easily,
 * and the inverse cdf computed once and cached.
 * If T2 is integral and f is "slow" in places, inv_f might have "holes".
 * For example, if f is erf, around 0 and 1, it's almost "flat", which will result
 * in inv_f being almost "vertical", which will result in missed indices
 * if T2 is integral.
 */
template<typename T1, typename T2>
inline void inverse_f(const std::vector<T1>& f, std::vector<T2>& inv_f) {
  for (size_t i = 1; i != f.size(); ++i)
    assert(f[i] != f[i - 1]);

  bool increasing = f[1] - f[0] > 0;

  if (increasing)
    for (size_t i = 1; i != f.size(); ++i)
      assert(f[i] > f[i - 1]);
  else
    for (size_t i = 1; i != f.size(); ++i)
      assert(f[i] < f[i - 1]);

  size_t n = inv_f.size();
  double f_min = min_value(f);
  double f_max = max_value(f);
  double y_step = (f_max - f_min) / n;
  size_t i = 0;

  for (double y = f_min; y <= f_max; y += y_step, ++i) {

    size_t idx = 0;
    if (increasing)
      while (idx + 1 < f.size() && f[idx + 1] < y)
        ++idx;
    else
      while (idx + 1 < f.size() && y < f[idx + 1])
        ++idx;

    double s = 1.0 / (double) (f[idx + 1] - f[idx]);
    double x = idx + s * (y - f[idx]);

    inv_f[i] = (T2) x;
  }
}

//--------------------------------------------------------------------------------
/**
 * Sample from a given cdf. From the ** CDF **, not the pdf!!
 */
template<typename T1, typename T2>
inline void sample(const std::vector<T1>& cdf, std::vector<T2>& samples) {
  {
    assert(!cdf.empty());
    assert(!samples.empty());
  }

  size_t n = samples.size();
  double M = cdf[cdf.size() - 1] / double(RAND_MAX);

  for (size_t i = 0; i != n; ++i) {
    double s = cdf[0], p = M * double(rand());
    size_t k = 0;
    while (s < p && k < cdf.size())
      s = cdf[++k];
    samples[i] = (T2) k;
  }
}

//--------------------------------------------------------------------------------
template<typename I>
inline void histogram(const std::vector<I>& samples, std::vector<I>& hist) {
  assert(!samples.empty());

  set_to_zero(hist);
  const I *i = &samples[0], *i_end = i + samples.size();

  for (; i != i_end; ++i) {
    if (hist.size() <= *i)
      hist.resize((*i) + 1, 0);
    ++hist[*i];
  }
}

//--------------------------------------------------------------------------------
/**
 * Generates a matrix of random (index,value) pairs from [0..ncols], with
 * nnzpr numbers per row, and n columns [generates a constant sparse matrix
 * with constant number of non-zeros per row]. This uses a 2D Gaussian distribution
 * for the on bits of each coincidence. That is, each coincidence is seen as a
 * folded 2D array, and a 2D Gaussian is used to distribute the on bits of each
 * coincidence.
 *
 * Each row is seen as an image of size (ncols / rf_x) by rf_x.
 * 'sigma' is the parameter of the Gaussian, which is centered at the center of
 * the image. Uses a symmetric Gaussian, specified only by the location of its
 * max and a single sigma parameter (no Sigma matrix). We use the symmetry of the
 * 2d gaussian to simplify computations. The conditional distribution obtained
 * from the 2d gaussian by fixing y is again a gaussian, with parameters than can
 * be easily deduced from the original 2d gaussian.
 */
template<typename T1, typename T2>
inline void
gaussian_2d_pair_sample(size_t nrows, size_t ncols, size_t nnzpr, size_t rf_x,
                        T2 sigma,
                        std::vector<std::pair<T1, T2> >& a,
                        const T2& init_nz_val,
                        int seed = -1,
                        bool sorted = true) {
  {
    assert(ncols % rf_x == 0);
    assert(nnzpr <= ncols);
    assert(0 < sigma);
  }

  a.resize(nrows * nnzpr);

  size_t rf_y = ncols / rf_x;
  T2 c_x = float(rf_x - 1.0) / 2.0, c_y = float(rf_y - 1.0) / 2.0;
  Gaussian2D<float> sg2d(c_x, c_y, sigma * sigma, 0, 0, sigma * sigma);
  std::vector<float> z(ncols);

  // Renormalize because we've lost some mass
  // with a compact domain of definition.
  float s = 0;
  for (size_t j = 0; j != ncols; ++j)
    s += z[j] = sg2d(j / rf_y, j % rf_y);

  for (size_t j = 0; j != ncols; ++j)
    z[j] /= s;
  // Test with uniform in case of bug or bias
  //z[j] = 1.0f / (float)(rf_x * rf_y);

  //std::vector<int> counts(ncols, 0);

  // TODO: argsort z so that the bigger bins come first, and it's faster
  // to draw samples in the area where the pdf is higher
  for (size_t i = 0; i != nrows; ++i) {

    std::set<size_t> b;

    while (b.size() < nnzpr) {
      T2 s = z[0], p = T2(rand()) / T2(RAND_MAX);
      size_t k = 0;
      while (s < p && k < ncols - 1)
        s += z[++k];
      //++counts[k];
      b.insert(k);
    }

    size_t offset = i * nnzpr;
    std::set<size_t>::const_iterator it = b.begin();
    for (size_t j = 0; j != nnzpr; ++j, ++it)
      a[offset + j] = std::make_pair<T1, T2>(*it, init_nz_val);
  }

  /*
    for (size_t i = 0; i != counts.size(); ++i)
    std::cout << counts[i] << " ";
    std::cout << std::endl;
  */
}

//--------------------------------------------------------------------------------
// Computes expectation from array that represents histogram of distribution.
// The categories of the histogram are implicitly 0,1,2,... hist.size()-1.
template<typename T>
inline double dist_expected(const std::vector<T>& hist) {
  { // Pre-conditions
#ifdef DEBUG
    assert(! hist.empty());
    double s = 0.0;
    for (size_t i = 0; i != hist.size(); ++i) {
      assert(0 <= hist[i]);
      s += hist[i];
    }
    assert(0 < s);
#endif
  } // End pre-conditions

  double sum = 0.0, mu = 0.0;

  for (size_t i = 0; i != hist.size(); ++i) {
    sum += hist[i];
    mu += i * hist[i];
  }

  mu /= sum;

  return mu;
}

//--------------------------------------------------------------------------------
/**
 * Finds percentiles.
 */
template<typename It1, typename It2>
inline void percentiles(size_t n_percentiles,
                        It1 begin1, It1 end1,
                        It2 begin2, It2 end2,
                        bool alreadyNormalized = false) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typedef typename std::iterator_traits<It1>::value_type value_type;
  typedef typename std::iterator_traits<It2>::value_type size_type;

  value_type n = (value_type) (alreadyNormalized ? 1.0f : 0.0f);

  if (!alreadyNormalized)
    for (It1 it = begin1; it != end1; ++it)
      n += *it;

  value_type increment = n / value_type(n_percentiles);
  value_type sum = (value_type) 0.0f;
  size_type p = (size_type) 0;

  for (value_type v = increment; v < n; v += increment) {
    for (; sum < v; ++p)
      sum += *begin1++;
    *begin2++ = p;
  }
}

//--------------------------------------------------------------------------------
template<typename C1, typename C2>
inline void percentiles(size_t n_percentiles, const C1& pmf, C2& pcts) {
  percentiles(n_percentiles, pmf.begin(), pmf.end(), pcts.begin());
}

//--------------------------------------------------------------------------------
template<typename T>
inline std::pair<double, double> mean_stddev(const std::vector<T>& x) {
  double s = sum(x);
  double avg = s / x.size();
  double ee = 0.0;
  for (size_t i = 0; i != x.size(); ++i)
    ee += (x[i] - avg) * (x[i] - avg);
  return std::make_pair(avg, sqrt(ee / (x.size() - 1)));
}

//--------------------------------------------------------------------------------
// PERMUTATIONS
//--------------------------------------------------------------------------------
// Based on Lehmer, factorial numbering, slow
/*
  template <typename T>
  struct permutation_iterator
  {
  size_t i, n;
  vector<size_t> indices;
  vector<T> p;

  inline permutation_iterator(const vector<T>& _p)
  : i(0), n(_p.size()), indices(n, 0), p(_p)
  {}

  inline permutation_iterator(const permutation_iterator& o)
  : i(o.i), n(o.n), indices(o.indices), p(o.p)
  {}

  inline permutation_iterator& operator=(const permutation_iterator& o)
  {
  if (&o != this) {
  i = o.i;
  n = o.n;
  indices = o.indices;
  p = o.p;
  }
  return *this;
  }

  inline vector<T> operator*()
  {
  vector<T> r;
  list<T> tmp;
  for (size_t j = 0; j != n; ++j)
  tmp.push_back(p[j]);
  for (int j = (int)n-1; 0 <= j; --j) {
  typename list<T>::iterator it = tmp.begin();
  for (size_t k = 0; k != indices[j]; ++k)
  ++it;
  r.push_back(*it);
  tmp.erase(it);
  }
  return r;
  }

  inline bool ok() const
  {
  return i < n;
  }

  inline void operator++()
  {
  do {
  if (indices[i] < i) {
  ++ indices[i];
  i = 0;
  } else {
  indices[i] = 0;
  ++ i;
  }
  } while (i != 0 && i < n);
  }
  };
*/

//--------------------------------------------------------------------------------
// Gray code for permutations.
// Algorithm following B.R.Heap "Permutations by interchanges" (1963)
// Optimized implementation, very fast.
template<typename T>
struct permutation_iterator {
  using ulong = unsigned long;
  ulong d_[32];
  ulong n_;   // permutations of n elements
  ulong ct_;  // count 5,4,3,2,1,(0); nonzero ==> easy cases
  std::vector<T> x;

  inline permutation_iterator(const std::vector<T>& _x)  // must have n>=3
  : n_(_x.size()), ct_(5), x(_x) {
    // d[0] and d[1] are unused
    ulong s = n_ < 3 ? 3 : n_;
    for (ulong k = 0; k < s; ++k)
      d_[k] = 0;
    d_[s - 1] = -1UL;  // sentinel
  }

  inline std::vector<T>& operator*() {
    return x;
  }

  // Return index of last element with reversal.
  // Return n with last permutation.
  inline bool next() {
    if (ct_ != 0) { // easy cases

      --ct_;
      ulong sw1 = 1 + (ct_ & 1);  // == 1,2,1,2,1
      std::swap(x[sw1], x[0]);
      return sw1 != n_;

    } else {

      ct_ = 5;  // reset counter

      // increment mixed radix number:
      ulong j = 2;
      while (d_[j] == j + 1)
        d_[j++] = 0;  // can touch sentinel

      // j == n-1 for last permutation:
      if (j == n_ - 1)
        return false;

      ulong k = j + 1;
      ulong xx = k & 1 ? d_[j] : 0;
      std::swap(x[k], x[xx]);

      ++d_[j];

      return k != n_;
    }
  }
};

//--------------------------------------------------------------------------------
// Watch out if not starting at the first one in the lexicographical order!
template<typename T>
inline bool next_permutation_lexi(T& x) {
  return std::next_permutation(x.begin(), x.end());
}

//--------------------------------------------------------------------------------
// RAND
//--------------------------------------------------------------------------------
template<typename T>
inline void random_shuffle(std::vector<T>& x) {
  std::random_shuffle(x.begin(), x.end());
}

//--------------------------------------------------------------------------------
template<typename It, typename RNG>
inline void rand_range(It begin, It end,
                       const typename std::iterator_traits<It>::value_type& min_,
                       const typename std::iterator_traits<It>::value_type& max_,
                       RNG& rng) {
  {
    assert(begin <= end);
    assert(min_ < max_);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  double range = double(max_ - min_) / double(rng.max() - rng.min());
  for (; begin != end; ++begin)
    *begin = value_type(double(rng()) * range + min_);
}

//--------------------------------------------------------------------------------
/**
 * Initializes a range with the normal distribution.
 *
 * @param begin
 * @param end
 * @param mean
 * @param stddev
 */
template<typename It>
inline void normal_range(It begin, It end,
                         const typename std::iterator_traits<It>::value_type& mean,
                         const typename std::iterator_traits<It>::value_type& stddev) {
  {
    assert(begin <= end);
  }

  //typedef typename std::iterator_traits<It>::value_type value_type;
  // implement numerical recipes' method
}

//--------------------------------------------------------------------------------
template<typename It, typename RNG>
inline void rand_range_01(It begin, It end, double pct, RNG& rng) {
  {
    assert(begin <= end);
    assert(0 <= pct && pct < 1);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  for (; begin != end; ++begin)
    *begin = (value_type) (double(rng()) / double(rng.max() - rng.min()) > pct);
}

//--------------------------------------------------------------------------------
template<typename T, typename RNG>
inline void rand_range_01(T& a, double pct, RNG& rng) {
  rand_range_01(a.begin(), a.end(), pct, rng);
}

//--------------------------------------------------------------------------------
/**
 * Initializes a range with a ramp function.
 *
 * @param begin
 * @param end
 * @param start the start value of the ramp
 * @param step the step of the ramp
 */
template<typename It, typename T>
inline void ramp_range(It begin, It end, T start = 0, T step = 1) {
  {
    assert(begin <= end);
  }

  for (; begin != end; ++begin, start += step)
    *begin = start;
}

//--------------------------------------------------------------------------------
/**
 * Initializes a range with a ramp function.
 *
 * @param a the container
 * @param start the start value of the ramp
 * @param step the step of the ramp
 */
template<typename T>
inline void ramp_range(T& a,
                       typename T::value_type start = 0,
                       typename T::value_type step = 1) {
  ramp_range(a.begin(), a.end(), start);
}

//--------------------------------------------------------------------------------
/**
 * Fills a range with values taken randomly from another range.
 *
 * @param begin
 * @param end
 * @param enum_begin the values to draw from
 * @param enum_end the values to draw from
 * @param replace whether to draw with or without replacements
 */
template<typename It1, typename It2, typename RNG>
inline void rand_enum_range(It1 begin, It1 end, It2 enum_begin, It2 enum_end,
                            bool replace, RNG& rng) {
  {
    assert(begin <= end);
    assert(enum_begin <= enum_end);
  }

  typedef typename std::iterator_traits<It1>::value_type value_type;

  size_t n = (size_t) (enum_end - enum_begin);

  if (replace) {

    for (; begin != end; ++begin)
      *begin = (value_type) *(enum_begin + rng() % n);

  } else {

    std::vector<size_t> ind(n);
    ramp_range(ind);

    for (; begin != end; ++begin) {
      size_t p = rng() % ind.size();
      *begin = (value_type) *(enum_begin + p);
      remove(p, ind);
    }
  }
}

//--------------------------------------------------------------------------------
/**
 * Fills a range with a random permutation of a ramp.
 *
 * @param begin
 * @param end
 * @param start the start value of the ramp
 * @param step the step value of the ramp
 */
template<typename It, typename RNG>
inline void
random_perm_interval(It begin, It end,
                     typename std::iterator_traits<It>::value_type start,
                     typename std::iterator_traits<It>::value_type step,
                     RNG& rng) {
  {
    assert(begin <= end);
  }

  ramp_range(begin, end, start, step);
  std::random_shuffle(begin, end, rng);
}

//--------------------------------------------------------------------------------
template<typename T>
inline void random_sample(size_t n, std::vector<T>& a, bool sorted = true) {
  assert(0 < a.size());

  std::vector<size_t> x(n);
  for (size_t i = 0; i != n; ++i)
    x[i] = i;
  std::random_shuffle(x.begin(), x.end());
  std::copy(x.begin(), x.begin() + a.size(), a.begin());
  if (sorted)
    std::sort(a.begin(), a.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline void random_sample(const std::set<T>& a, std::vector<T>& b, bool sorted = true) {
  assert(0 < b.size());

  std::vector<T> aa(a.begin(), a.end());
  std::random_shuffle(aa.begin(), aa.end());
  std::copy(aa.begin(), aa.begin() + b.size(), b.begin());
  if (sorted)
    std::sort(b.begin(), b.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline void random_binary(float proba, std::vector<T>& x) {
  size_t threshold = (size_t) (proba * 65535);
  for (size_t i = 0; i != x.size(); ++i)
    x[i] = ((size_t) rand() % 65535 < threshold) ? 1 : 0;
}

////--------------------------------------------------------------------------------
///**
// * Generates a matrix of random indices samples from [0..ncols], with
// * nnzpc numbers per row, and n columns [generates the indices for a sparse binary
// * matrix with constant number of non-zeros per row].
// */
//template <typename T>
//inline void
//random_sample(size_t nrows, size_t ncols, size_t nnzpr, std::vector<T>& a,
//              bool sorted =true)
//{
//  assert(0 < a.size());
//  assert(a.size() >= nrows * nnzpr);
//
//  Random rng; // uses srand/rand
//
//  std::vector<size_t> x(ncols);
//  for (size_t i = 0; i != ncols; ++i)
//    x[i] = i;
//  for (size_t i = 0; i != nrows; ++i) {
//    std::random_shuffle(x.begin(), x.end(), rng);
//    if (sorted)
//      std::sort(x.begin(), x.begin() + nnzpr);
//    std::copy(x.begin(), x.begin() + nnzpr, a.begin() + i * nnzpr);
//  }
//}

////--------------------------------------------------------------------------------
///**
// * Generates a matrix of random (index,value) pairs from [0..ncols], with
// * nnzpc numbers per row, and n columns [generates a constant sparse matrix
// * with constant number of non-zeros per row].
// */
//template <typename T1, typename T2>
//inline void
//random_pair_sample(size_t nrows, size_t ncols, size_t nnzpr,
//                   std::vector<std::pair<T1, T2> >& a,
//                   const T2& init_nz_val,
//                   bool sorted =true)
//{
//  assert(0 < a.size());
//
//  a.resize(nrows * nnzpr);
//
//  Random rng; // uses srand/rand
//
//  std::vector<size_t> x(ncols);
//  for (size_t i = 0; i != ncols; ++i)
//    x[i] = i;
//  for (size_t i = 0; i != nrows; ++i) {
//    std::random_shuffle(x.begin(), x.end(), rng);
//    if (sorted)
//      std::sort(x.begin(), x.begin() + nnzpr);
//    size_t offset = i*nnzpr;
//    for (size_t j = 0; j != nnzpr; ++j)
//      a[offset + j] = std::make_pair<T1,T2>(x[j], init_nz_val);
//  }
//}

//--------------------------------------------------------------------------------
/**
 * Generates a matrix of random (index,value) pairs from [0..ncols], with
 * nnzpc numbers per row, and n columns [generates a random sparse matrix, that
 * if a sparse matrix of random values, with a fixed number of non-zeros per row].
 */
template<typename T1, typename T2, typename RNG>
inline void
random_pair_sample_rng(size_t nrows, size_t ncols, size_t nnzpr,
                       std::vector<std::pair<T1, T2> >& a,
                       RNG& rng,
                       bool sorted = true) {
  assert(0 < a.size());

  a.resize(nrows * nnzpr);

  std::vector<size_t> x(ncols);
  for (size_t i = 0; i != ncols; ++i)
    x[i] = i;
  for (size_t i = 0; i != nrows; ++i) {
    std::random_shuffle(x.begin(), x.end(), rng);
    if (sorted)
      std::sort(x.begin(), x.begin() + nnzpr);
    size_t offset = i * nnzpr;
    for (size_t j = 0; j != nnzpr; ++j)
      a[offset + j] = std::make_pair<T1, T2>(x[j], rng());
  }
}

//--------------------------------------------------------------------------------
// generate
//--------------------------------------------------------------------------------
/**
 * Initializes a range by calling gen() repetitively.
 *
 * @param c the container to initialize
 * @param gen the generator functor
 */
template<typename Container, typename Generator>
inline void generate(Container& c, Generator gen) {
  typename Container::iterator i = c.begin(), e = c.end();

  for (; i != e; ++i)
    *i = gen();
}

//--------------------------------------------------------------------------------
// concatenate
//--------------------------------------------------------------------------------
/**
 * Concatenates multiple sub-ranges of a source range into a single range.
 *
 * @param x_begin the beginning of the source range
 * @param seg_begin the beginning of the ranges that describe the
 *  sub-ranges (start, size)
 * @param seg_end one past the end of the ranges that describe the
 *  sub-ranges
 * @param y_begin the beginning of the concatenated range
 */
template<typename InIt1, typename InIt2, typename OutIt>
inline void
concatenate(InIt1 x_begin, InIt2 seg_begin, InIt2 seg_end, OutIt y_begin) {
  {
    assert(seg_begin <= seg_end);
  }

  for (; seg_begin != seg_end; ++seg_begin) {
    InIt1 begin = x_begin + seg_begin->first;
    InIt1 end = begin + seg_begin->second;
    std::copy(begin, end, y_begin);
    y_begin += seg_begin->second;
  }
}

//--------------------------------------------------------------------------------
// Clip, threshold, binarize
//--------------------------------------------------------------------------------
/**
 * Clip the values in a range to be between min (included) and max (included):
 * any value less than min becomes min and any value greater than max becomes
 * max.
 *
 * @param begin
 * @param end
 * @param _min the minimum value
 * @param _max the maximum value
 */
template<typename It>
inline void clip(It begin, It end,
                 const typename std::iterator_traits<It>::value_type& _min,
                 const typename std::iterator_traits<It>::value_type& _max) {
  {
    assert(begin <= end);
  }

  while (begin != end) {
    typename std::iterator_traits<It>::value_type val = *begin;
    if (val > _max)
      *begin = _max;
    else if (val < _min)
      *begin = _min;
    ++begin;
  }
}

//--------------------------------------------------------------------------------
/**
 * Clip the values in a container to be between min (included) and max (included):
 * any value less than min becomes min and any value greater than max becomes
 * max.
 *
 * @param a the container
 * @param _min
 * @param _max
 */
template<typename T>
inline void clip(T& a,
                 const typename T::value_type& _min,
                 const typename T::value_type& _max) {
  clip(a.begin(), a.end(), _min, _max);
}

//--------------------------------------------------------------------------------
/**
 * Threshold a range and puts the values that were not eliminated into
 * another (sparse) range (index, value).
 *
 * @param begin
 * @param end
 * @param ind the beginning of the sparse indices
 * @param nz the beginning of the sparse values
 * @param th the threshold to use
 */
template<typename InIter, typename OutIter1, typename OutIter2>
inline size_t threshold(InIter begin, InIter end,
                        OutIter1 ind, OutIter2 nz,
                        const typename std::iterator_traits<InIter>::value_type& th,
                        bool above = true) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<InIter>::value_type value_type;
  typedef size_t size_type;

  size_type n = 0;

  if (above) {

    for (InIter it = begin; it != end; ++it) {
      value_type val = (value_type) *it;
      if (val >= th) {
        *ind = (size_type) (it - begin);
        *nz = val;
        ++ind;
        ++nz;
        ++n;
      }
    }

  } else {

    for (InIter it = begin; it != end; ++it) {
      value_type val = (value_type) *it;
      if (val < th) {
        *ind = (size_type) (it - begin);
        *nz = val;
        ++ind;
        ++nz;
        ++n;
      }
    }
  }

  return n;
}

//--------------------------------------------------------------------------------
/**
 * Given a threshold and a dense vector x, returns another dense vectors with 1's
 * where the value of x is > threshold, and 0 elsewhere. Also returns the count
 * of 1's.
 */
template<typename InputIterator, typename OutputIterator>
inline unsigned int
binarize_with_threshold(float threshold,
                        InputIterator x, InputIterator x_end,
                        OutputIterator y, OutputIterator y_end) {
  {
    assert(x_end - x == y_end - y);
  }

  unsigned int count = 0;

  for (; x != x_end; ++x, ++y)
    if (*x > threshold) {
      *y = 1;
      ++count;
    } else
      *y = 0;

  return count;
}

//--------------------------------------------------------------------------------
// INDICATORS
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
/**
 * Given a dense 2D array of 0 and 1, return a vector that has as many rows as x
 * a 1 wherever x as a non-zero row, and a 0 elsewhere. I.e. the result is the
 * indicator of non-zero rows. Gets fast by not scanning a row more than is
 * necessary, i.e. stops as soon as a 1 is found on the row.
 */
template<typename InputIterator, typename OutputIterator>
inline void
nonZeroRowsIndicator_01(unsigned int nrows, unsigned int ncols,
                        InputIterator x, InputIterator x_end,
                        OutputIterator y, OutputIterator y_end) {
  {
    assert(0 < nrows);
    assert(0 < ncols);
    assert((unsigned int) (x_end - x) == nrows * ncols);
    assert((unsigned int) (y_end - y) == nrows);
#ifdef ASSERTIONS_ON
    for (unsigned int i = 0; i != nrows * ncols; ++i)
      assert(x[i] == 0 || x[i] == 1);
#endif
  }

  for (unsigned int r = 0; r != nrows; ++r, ++y) {

    InputIterator it = x + r * ncols, it_end = it + ncols;
    unsigned int found = 0;

    while (it != it_end && found == 0)
      found = (unsigned int) (*it++);

    *y = found;
  }
}

//--------------------------------------------------------------------------------
/**
 * Given a dense 2D array of 0 and 1, return the number of rows that have
 * at least one non-zero. Gets fast by not scanning a row more than is
 * necessary, i.e. stops as soon as a 1 is found on the row.
 */
template<typename InputIterator>
inline unsigned int
nNonZeroRows_01(unsigned int nrows, unsigned int ncols,
                InputIterator x, InputIterator x_end) {
  {
    assert(0 < nrows);
    assert(0 < ncols);
    assert((unsigned int) (x_end - x) == nrows * ncols);
#ifdef ASSERTIONS_ON
    for (unsigned int i = 0; i != nrows * ncols; ++i)
      assert(x[i] == 0 || x[i] == 1);
#endif
  }

  unsigned int count = 0;

  for (unsigned int r = 0; r != nrows; ++r) {

    InputIterator it = x + r * ncols, it_end = it + ncols;
    unsigned int found = 0;

    while (it != it_end && found == 0)
      found = (unsigned int) (*it++);

    count += found;
  }

  return count;
}

//--------------------------------------------------------------------------------
/**
 * Given a dense 2D array of 0 and 1 x, return a vector that has as many cols as x
 * a 1 wherever x as a non-zero col, and a 0 elsewhere. I.e. the result is the
 * indicator of non-zero cols. Gets fast by not scanning a row more than is
 * necessary, i.e. stops as soon as a 1 is found on the col.
 */
template<typename InputIterator, typename OutputIterator>
inline void
nonZeroColsIndicator_01(unsigned int nrows, unsigned int ncols,
                        InputIterator x, InputIterator x_end,
                        OutputIterator y, OutputIterator y_end) {
  {
    assert(0 < nrows);
    assert(0 < ncols);
    assert((unsigned int) (x_end - x) == nrows * ncols);
    assert((unsigned int) (y_end - y) == ncols);
#ifdef ASSERTIONS_ON
    for (unsigned int i = 0; i != nrows * ncols; ++i)
      assert(x[i] == 0 || x[i] == 1);
#endif
  }

  unsigned int N = nrows * ncols;

  for (unsigned int c = 0; c != ncols; ++c, ++y) {

    InputIterator it = x + c, it_end = it + N;
    unsigned int found = 0;

    while (it != it_end && found == 0) {
      found = (unsigned int) (*it);
      it += ncols;
    }

    *y = found;
  }
}

//--------------------------------------------------------------------------------
/**
 * Given a dense 2D array of 0 and 1, return the number of columns that have
 * at least one non-zero.
 * Gets fast by not scanning a col more than is necessary, i.e. stops as soon as
 * a 1 is found on the col.
 */
template<typename InputIterator>
inline unsigned int
nNonZeroCols_01(unsigned int nrows, unsigned int ncols,
                InputIterator x, InputIterator x_end) {
  {
    assert(0 < nrows);
    assert(0 < ncols);
    assert((unsigned int) (x_end - x) == nrows * ncols);
#ifdef ASSERTIONS_ON
    for (unsigned int i = 0; i != nrows * ncols; ++i)
      assert(x[i] == 0 || x[i] == 1);
#endif
  }

  unsigned int count = 0;
  unsigned int N = nrows * ncols;

  for (unsigned int c = 0; c != ncols; ++c) {

    InputIterator it = x + c, it_end = it + N;
    unsigned int found = 0;

    while (it != it_end && found == 0) {
      found = (unsigned int) (*it);
      it += ncols;
    }

    count += found;
  }

  return count;
}

//--------------------------------------------------------------------------------
// MASK
//--------------------------------------------------------------------------------
/**
 * Mask an array.
 */
template<typename InIter>
inline void mask(InIter begin, InIter end, InIter zone_begin, InIter zone_end,
                 const typename std::iterator_traits<InIter>::value_type& v = 0,
                 bool maskOutside = true) {
  { // Pre-conditions
    assert(begin <= end);
    assert(zone_begin <= zone_end);
    assert(begin <= zone_begin && zone_end <= end);
  } // End pre-conditions

  //typedef typename std::iterator_traits<InIter>::value_type value_type;

  if (maskOutside) {
    std::fill(begin, zone_begin, v);
    std::fill(zone_end, end, v);
  } else {
    std::fill(zone_begin, zone_end, v);
  }
}

//--------------------------------------------------------------------------------
template<typename value_type>
inline void mask(std::vector<value_type>& x,
                 typename std::vector<value_type>::size_type zone_begin,
                 typename std::vector<value_type>::size_type zone_end,
                 const value_type& v = 0,
                 bool maskOutside = true) {
  { // Pre-conditions
    assert(0 <= zone_begin && zone_begin <= zone_end && zone_end <= x.size());
  } // End pre-conditions

  mask(x.begin(), x.end(), x.begin() + zone_begin, x.begin() + zone_end,
       v, maskOutside);
}

//--------------------------------------------------------------------------------
template<typename value_type1, typename value_type2>
inline void mask(std::vector<value_type1>& x, const std::vector<value_type2>& mask,
                 bool multiplyYesNo = false, value_type2 eps = (value_type2) 1e-6) {
  { // Pre-conditions
    assert(x.size() == mask.size());
  } // End pre-conditions

  typedef typename std::vector<value_type1>::size_type size_type;

  if (multiplyYesNo) {
    for (size_type i = 0; i != x.size(); ++i)
      if (!nearlyZero(mask[i]), eps)
        x[i] *= (value_type1) mask[i];
      else
        x[i] = (value_type1) 0;

  } else {
    for (size_type i = 0; i != x.size(); ++i)
      if (nearlyZero(mask[i], eps))
        x[i] = (value_type1) 0;
  }
}

//--------------------------------------------------------------------------------
// NORMS
//--------------------------------------------------------------------------------
/**
 * A class that provides init and operator(), to be used in distance
 * computations when using the Hamming (L0) norm.
 */
template<typename T>
struct Lp0 {
  typedef T value_type;

  inline value_type operator()(value_type& a, value_type b) const {
    value_type inc = value_type(b < -1e-6 || b > 1e-6);
    a += inc;
    return inc;
  }

  inline value_type root(value_type x) const { return x; }
};

//--------------------------------------------------------------------------------
/**
 * A class that provides init and operator(), to be used in distance
 * computations when using the Manhattan (L1) norm.
 */
template<typename T>
struct Lp1 {
  typedef T value_type;

  inline value_type operator()(value_type& a, value_type b) const {
    value_type inc = fabs(b); //b > 0.0 ? b : -b;
    a += inc;
    return inc;
  }

  inline value_type root(value_type x) const { return x; }
};

//--------------------------------------------------------------------------------
/**
 * A class that provides square and square root methods, to be
 * used in distance computations when using L2 norm.
 */
template<typename T>
struct Lp2 {
  typedef T value_type;

  Sqrt<value_type> s;

  inline value_type operator()(value_type& a, value_type b) const {
    value_type inc = b * b;
    a += inc;
    return inc;
  }

  inline value_type root(value_type x) const {
    return s(x);
  }
};

//--------------------------------------------------------------------------------
/**
 * A class that provides power p and root p methods, to be
 * used in distance computations using Lp norm.
 */
template<typename T>
struct Lp {
  typedef T value_type;

  Pow<value_type> pf;

  Lp(value_type p_)
  : p(p_), inv_p((value_type) 1.0) {
    // We allow only positive values of p for now, as this
    // keeps the root function monotonically increasing, which
    // results in further speed-ups.
    assert(p_ > (value_type) 0.0);

    inv_p = (value_type) 1.0 / p;
  }

  value_type p, inv_p;

  inline value_type operator()(value_type& a, value_type b) const {
    value_type inc = pf(b > 0.0 ? b : -b, p);
    a += inc;
    return inc;
  }

  inline value_type root(value_type x) const {
    // skipping abs, because we know we've been adding positive
    // numbers when root is called when computing a norm
    return pf(x, inv_p);
  }
};

//--------------------------------------------------------------------------------
/**
 * A class that provides power p and root p methods, to be
 * used in distance computations using LpMax norm.
 */
template<typename T>
struct LpMax {
  typedef T value_type;

  Max<value_type> m;

  inline value_type operator()(value_type& a, value_type b) const {
    value_type inc = m(a, b > 0 ? b : -b);
    a = inc;
    return inc;
  }

  inline value_type root(value_type x) const { return x; }
};

//--------------------------------------------------------------------------------
/**
 * Hamming norm.
 */
template<typename It>
inline typename std::iterator_traits<It>::value_type
l0_norm(It begin, It end, bool = true) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  value_type n = (value_type) 0;
  Lp0 <value_type> lp0;

  for (; begin != end; ++begin)
    lp0(n, *begin);

  return n;
}

//--------------------------------------------------------------------------------
/**
 * Hamming norm on a container
 */
template<typename T>
inline typename T::value_type l0_norm(const T& a, bool = true) {
  return l0_norm(a.begin(), a.end());
}

//--------------------------------------------------------------------------------
/**
 * Manhattan norm.
 */
template<typename It>
inline typename std::iterator_traits<It>::value_type
l1_norm(It begin, It end, bool = true) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  value_type n = (value_type) 0;
  Lp1 <value_type> lp1;

  for (; begin != end; ++begin)
    lp1(n, *begin);

  return n;
}

//--------------------------------------------------------------------------------
/**
 * Manhattan norm on a container.
 */
template<typename T>
inline typename T::value_type l1_norm(const T& a, bool = true) {
  return l1_norm(a.begin(), a.end());
}

//--------------------------------------------------------------------------------
/**
 * Euclidean norm.
 */

//--------------------------------------------------------------------------------
#ifdef PLATFORM_darwin86
inline void sum_of_squares(float* begin, int n, float* s)
{
  vDSP_svesq(begin, 1, s, n);
}

//--------------------------------------------------------------------------------
inline void sum_of_squares(double* begin, int n, double* s)
{
  vDSP_svesqD(begin, 1, s, n);
}
#endif

//--------------------------------------------------------------------------------
template<typename It>
inline typename std::iterator_traits<It>::value_type
l2_norm(It begin, It end, bool take_root = true) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;
  value_type n = (value_type) 0;

  Lp2 <value_type> lp2;

#ifdef PLATFORM_darwin86 // 10X faster

  // &*begin won't work on platforms where the iterators are not pointers
  // (win32)
  sum_of_squares(&*begin, (end - begin), &n);

#else

  for (; begin != end; ++begin)
    lp2(n, *begin);

#endif

  if (take_root)
    n = lp2.root(n);

  return n;
}

//--------------------------------------------------------------------------------
/**
 * Euclidean norm on a container.
 */
template<typename T>
inline typename T::value_type l2_norm(const T& a, bool take_root = true) {
  return l2_norm(a.begin(), a.end(), take_root);
}

//--------------------------------------------------------------------------------
/**
 * p-norm.
 */
template<typename It>
inline typename std::iterator_traits<It>::value_type
lp_norm(typename std::iterator_traits<It>::value_type p,
        It begin, It end, bool take_root = true) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  value_type n = (value_type) 0;
  Lp <value_type> lp(p);

  for (; begin != end; ++begin)
    lp(n, *begin);

  if (take_root)
    n = lp.root(n);

  return n;
}

//--------------------------------------------------------------------------------
/**
 * p-norm on a container.
 */
template<typename T>
inline typename T::value_type
lp_norm(typename T::value_type p, const T& a, bool take_root = true) {
  return lp_norm(p, a.begin(), a.end(), take_root);
}

//--------------------------------------------------------------------------------
/**
 * L inf / L max norm.
 */
template<typename It>
inline typename std::iterator_traits<It>::value_type
lmax_norm(It begin, It end, bool = true) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  value_type n = (value_type) 0;
  LpMax <value_type> lmax;

  for (; begin != end; ++begin)
    lmax(n, *begin);

  return n;
}

//--------------------------------------------------------------------------------
/**
 * L inf / L max norm on a container.
 */
template<typename T>
inline typename T::value_type lmax_norm(const T& a, bool = true) {
  return lmax_norm(a.begin(), a.end());
}

//--------------------------------------------------------------------------------
/**
 * Norm function.
 *
 * @param p the norm
 * @param begin
 * @param end
 * @param take_root whether to take the p-th root or not
 */
template<typename It>
inline typename std::iterator_traits<It>::value_type
norm(typename std::iterator_traits<It>::value_type p,
     It begin, It end, bool take_root = true) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  if (p == (value_type) 0)
    return l0_norm(begin, end);
  else if (p == (value_type) 1)
    return l1_norm(begin, end);
  else if (p == (value_type) 2)
    return l2_norm(begin, end, take_root);
  else if (p == std::numeric_limits<value_type>::max())
    return lmax_norm(begin, end);
  else
    return lp_norm(p, begin, end, take_root);
}

//--------------------------------------------------------------------------------
/**
 * Norm on a whole container.
 */
template<typename T>
inline typename T::value_type
norm(typename T::value_type p, const T& a, bool take_root = true) {
  return norm(p, a.begin(), a.end(), take_root);
}

//--------------------------------------------------------------------------------
/**
 * Normalize a range, according to the p norm, so that the values sum up to n.
 *
 * @param begin
 * @param end
 * @param p the norm
 * @param n the value of the sum of the elements after normalization
 */
template<typename It>
inline void
normalize(It begin, It end,
          const typename std::iterator_traits<It>::value_type& p = 1.0,
          const typename std::iterator_traits<It>::value_type& n = 1.0) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  value_type s = (value_type) 0;

  if (p == (value_type) 0)
    s = l0_norm(begin, end);
  else if (p == (value_type) 1)
    s = l1_norm(begin, end);
  else if (p == (value_type) 2)
    s = l2_norm(begin, end);
  else if (p == std::numeric_limits<value_type>::max())
    s = lmax_norm(begin, end);

  if (s != (value_type) 0)
    multiply_val(begin, end, n / s);
}

//--------------------------------------------------------------------------------
/**
 * Normalize a container, with p-th norm and so that values add up to n.
 */
template<typename T>
inline void normalize(T& a,
                      const typename T::value_type& p = 1.0,
                      const typename T::value_type& n = 1.0) {
  normalize(a.begin(), a.end(), p, n);
}

//--------------------------------------------------------------------------------
/**
 * Normalization according to LpMax: finds the max of the range,
 * and then divides all the values so that the max is n.
 * Makes it nicer to call normalize when using LpMax.
 */
template<typename It>
inline void
normalize_max(It begin, It end,
              const typename std::iterator_traits<It>::value_type& n = 1.0) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;

  normalize(begin, end, std::numeric_limits<value_type>::max(), n);
}

//--------------------------------------------------------------------------------
/**
 * Normalization according to LpMax.
 */
template<typename value_type>
inline void normalize_max(std::vector<value_type>& x, const value_type& n = 1.0) {
  normalize_max(x.begin(), x.end(), n);
}

//--------------------------------------------------------------------------------
// DISTANCES
//--------------------------------------------------------------------------------
/**
 * Returns the max of the absolute values of the differences.
 *
 * @param begin1
 * @param end1
 * @param begin2
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
max_abs_diff(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typename std::iterator_traits<It1>::value_type d(0), val(0);

  while (begin1 != end1) {
    val = *begin1 - *begin2;
    val = val > 0 ? val : -val;
    if (val > d)
      d = val;
    ++begin1;
    ++begin2;
  }

  return d;
}

//--------------------------------------------------------------------------------
/**
 * Returns the max of the absolute values of the differences.
 *
 * @param a first container
 * @param b second container
 */
template<typename T1, typename T2>
inline typename T1::value_type max_abs_diff(const T1& a, const T2& b) {
  return max_abs_diff(a.begin(), a.end(), b.begin(), b.end());
}

//--------------------------------------------------------------------------------
/**
 * Returns the Hamming distance of the two ranges.
 *
 * @param begin1
 * @param end1
 * @param begin2
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
hamming_distance(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typename std::iterator_traits<It1>::value_type d(0);

  while (begin1 != end1) {
    d += *begin1 != *begin2;
    ++begin1;
    ++begin2;
  }

  return d;
}

//--------------------------------------------------------------------------------
/**
 * Returns the Hamming distance of the two containers.
 *
 * @param a first container
 * @param b second container
 */
template<typename T1, typename T2>
inline typename T1::value_type
hamming_distance(const T1& a, const T2& b) {
  return hamming_distance(a.begin(), a.end(), b.begin(), b.end());
}

//--------------------------------------------------------------------------------
/**
 * [begin1, end1) and [begin2, end2) are index encodings of binary 0/1 ranges.
 */
template<typename It1, typename It2>
inline size_t
sparse_hamming_distance(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    // todo: check that ranges are valid sparse indices ranges
    // (increasing, no duplicate...)
    assert(begin1 <= end1);
    assert(begin2 <= end2);
  }

  typedef size_t size_type;

  size_type d = 0;

  while (begin1 != end1 && begin2 != end2) {
    if (*begin1 < *begin2) {
      ++d;
      ++begin1;
    } else if (*begin2 < *begin1) {
      ++d;
      ++begin2;
    } else {
      ++begin1;
      ++begin2;
    }
  }

  d += (size_type) (end1 - begin1);
  d += (size_type) (end2 - begin2);

  return d;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline typename T1::size_type
sparse_hamming_distance(const T1& a, const T2& b) {
  return sparse_hamming_distance(a.begin(), a.end(), b.begin(), b.end());
}

//--------------------------------------------------------------------------------
/**
 * Returns the Manhattan distance of the two ranges.
 *
 * @param begin1
 * @param end1
 * @param begin2
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
manhattan_distance(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typedef typename std::iterator_traits<It1>::value_type value_type;

  value_type d = (value_type) 0;
  Lp1 <value_type> lp1;

  for (; begin1 != end1; ++begin1, ++begin2)
    lp1(d, *begin1 - *begin2);

  return d;
}

//--------------------------------------------------------------------------------
/**
 * Returns the Manhattan distance of the two containers.
 *
 * @param a first container
 * @param b second container
 */
template<typename T1, typename T2>
inline typename T1::value_type
manhattan_distance(const T1& a, const T2& b) {
  return manhattan_distance(a.begin(), a.end(), b.begin(), b.end());
}

//--------------------------------------------------------------------------------
/**
 * Returns the Euclidean distance of the two ranges.
 *
 * @param begin1
 * @param end1
 * @param begin2
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
euclidean_distance(It1 begin1, It1 end1, It2 begin2, It2 end2, bool take_root = true) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typedef typename std::iterator_traits<It1>::value_type value_type;

  value_type d = (value_type) 0;
  Lp2 <value_type> lp2;

  for (; begin1 != end1; ++begin1, ++begin2)
    lp2(d, *begin1 - *begin2);

  if (take_root)
    d = lp2.root(d);

  return d;
}

//--------------------------------------------------------------------------------
/**
 * Returns the Euclidean distance of the two containers.
 *
 * @param a first container
 * @param b second container
 */
template<typename T1, typename T2>
inline typename T1::value_type
euclidean_distance(const T1& a, const T2& b, bool take_root = true) {
  return euclidean_distance(a.begin(), a.end(), b.begin(), b.end(), take_root);
}

//--------------------------------------------------------------------------------
/**
 * Returns the Lp distance of the two ranges.
 *
 * @param begin1
 * @param end1
 * @param begin2
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
lp_distance(typename std::iterator_traits<It1>::value_type p,
            It1 begin1, It1 end1, It2 begin2, It2 end2, bool take_root = true) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typedef typename std::iterator_traits<It1>::value_type value_type;

  value_type d = (value_type) 0;
  Lp <value_type> lp(p);

  for (; begin1 != end1; ++begin1, ++begin2)
    lp(d, *begin1 - *begin2);

  if (take_root)
    d = lp.root(d);

  return d;
}

//--------------------------------------------------------------------------------
/**
 * Returns the Lp distance of the two containers.
 *
 * @param a first container
 * @param b second container
 */
template<typename T1, typename T2>
inline typename T1::value_type
lp_distance(typename T1::value_type p,
            const T1& a, const T2& b, bool take_root = true) {
  return lp_distance(p, a.begin(), a.end(), b.begin(), b.end(), take_root);
}

//--------------------------------------------------------------------------------
/**
 * Returns the Lmax distance of the two ranges.
 *
 * @param begin1
 * @param end1
 * @param begin2
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
lmax_distance(It1 begin1, It1 end1, It2 begin2, It2 end2, bool = true) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typedef typename std::iterator_traits<It1>::value_type value_type;

  value_type d = (value_type) 0;
  LpMax <value_type> lmax;

  for (; begin1 != end1; ++begin1, ++begin2)
    lmax(d, *begin1 - *begin2);

  return d;
}

//--------------------------------------------------------------------------------
/**
 * Returns the Lmax distance of the two containers.
 *
 * @param a first container
 * @param b second container
 */
template<typename T1, typename T2>
inline typename T1::value_type
lmax_distance(const T1& a, const T2& b, bool = true) {
  return lmax_distance(a.begin(), a.end(), b.begin(), b.end());
}

//--------------------------------------------------------------------------------
/**
 * Returns the distance of the two ranges.
 *
 * @param begin1
 * @param end1
 * @param begin2
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
distance(typename std::iterator_traits<It1>::value_type p,
         It1 begin1, It1 end1, It2 begin2, It2 end2, bool take_root = true) {
  {
    assert(begin1 <= end1);
    assert(begin2 <= end2);
    assert(end1 - begin1 == end2 - begin2);
  }

  typedef typename std::iterator_traits<It1>::value_type value_type;

  if (p == (value_type) 0)
    return hamming_distance(begin1, end1, begin2);
  else if (p == (value_type) 1)
    return manhattan_distance(begin1, end1, begin2);
  else if (p == (value_type) 2)
    return euclidean_distance(begin1, end1, begin2, take_root);
  else if (p == std::numeric_limits<value_type>::max())
    return lmax_distance(begin1, end1, begin2);
  else
    return lp_distance(p, begin1, end1, begin2, take_root);
}

//--------------------------------------------------------------------------------
/**
 * Returns the distance of the two containers.
 *
 * @param a first container
 * @param b second container
 */
template<typename T1, typename T2>
inline typename T1::value_type
distance(typename T1::value_type p, const T1& a, const T2& b, bool take_root = true) {
  return distance(p, a.begin(), a.end(), b.begin(), b.end(), take_root);
}

//--------------------------------------------------------------------------------
// Counting
//--------------------------------------------------------------------------------
/**
 * Counts the elements which satisfy the passed predicate in the given range.
 */
template<typename C, typename Predicate>
inline size_t count_if(const C& c, Predicate pred) {
  return std::count_if(c.begin(), c.end(), pred);
}

//--------------------------------------------------------------------------------
/**
 * Counts the number of zeros in the given range.
 */
template<typename It>
inline size_t
count_zeros(It begin, It end,
            const typename std::iterator_traits<It>::value_type& eps = 1e-6) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;
  return std::count_if(begin, end, IsNearlyZero < DistanceToZero < value_type > > (eps));
}

//--------------------------------------------------------------------------------
/**
 * Counts the number of zeros in the container passed in.
 */
template<typename C>
inline size_t count_zeros(const C& c, const typename C::value_type& eps = 1e-6) {
  return count_zeros(c.begin, c.end(), eps);
}

//--------------------------------------------------------------------------------
/**
 * Count the number of ones in the given range.
 */
template<typename It>
inline size_t
count_ones(It begin, It end,
           const typename std::iterator_traits<It>::value_type& eps = 1e-6) {
  {
    assert(begin <= end);
  }

  typedef typename std::iterator_traits<It>::value_type value_type;
  return std::count_if(begin, end, IsNearlyZero < DistanceToOne < value_type > > (eps));
}

//--------------------------------------------------------------------------------
/**
 * Count the number of ones in the container passed in.
 */
template<typename C>
inline size_t count_ones(const C& c, const typename C::value_type& eps = 1e-6) {
  return count_ones(c.begin(), c.end(), eps);
}

//--------------------------------------------------------------------------------
/**
 * Counts the number of values greater than a given threshold in a given range.
 *
 * Asm SSE is many times faster than C++ (almost 10X), and C++ is 10X faster than
 * numpy (some_array > threshold).sum(). The asm code doesn't have branchs, which
 * is probably very good for the CPU front-end.
 *
 * This is not as general as a count_gt that would be parameterized on the type
 * of the elements in the range, and it requires passing in a Python arrays
 * that are .astype(float).
 *
 * Doesn't work for 64 bits platforms, doesn't work on win32.
 */
inline unsigned int
count_gt(float *begin, float *end, float threshold) {
  // Need this, because the asm syntax is not correct for win32,
  // we simply can't compile the code as is on win32.
#ifdef PLATFORM_darwin86

  // Need this, because even on darwin86, some older machines might
  // not have the right SSE instructions.
  if (SSE_LEVEL >= 3) {

    // Compute offsets into array [begin..end):
    // start is the first 4 bytes aligned address (to start movaps)
    // n0 is the number of floats before we reach start and can use parallel
    //  xmm operations
    // n1 is the number floats we can process in parallel with xmm
    // n2 is the number of "stragglers" what we will have to do one by one ( < 4)
    float count = 0;
    long x_addr = (long) begin; // 8 bytes on 64 bits platforms
    float* start = (x_addr % 16 == 0) ? begin : (float*) (16*(x_addr/16+1));
    int n0 = (int)(start - begin);
    int n1 = 4 * ((end - start) / 4);
    int n2 = (int)(end - start - n1);

    asm volatile(
                 // Prepare various xmm registers, storing the value of the
                 // threshold and the value 1: with xmm, we will operate on
                 // 4 floats at a time, so we replicate threshold 4 times in
                 // xmm1, and 4 times again in xmm2. The operations will be in
                 // parallel.
                 "subl $16, %%esp\n\t"            // allocate 4 floats on stack

                 "movl %%eax, (%%esp)\n\t"        // copy threshold to 4 locations
                 "movl %%eax, 4(%%esp)\n\t"       // on stack: we want threshold
                 "movl %%eax, 8(%%esp)\n\t"       // to be filling xmm1 and xmm2
                 "movl %%eax, 12(%%esp)\n\t"      // (operate on 4 floats at a time)
                 "movaps (%%esp), %%xmm1\n\t"     // move 4 thresholds into xmm1
                 "movaps %%xmm1, %%xmm2\n\t"      // copy 4 thresholds to xmm2

                 "movl $0x3f800000, (%%esp)\n\t"  // $0x3f800000 = (float) 1.0
                 "movl $0x3f800000, 4(%%esp)\n\t" // we want to have that constant
                 "movl $0x3f800000, 8(%%esp)\n\t" // 8 times, in xmm3 and xmm4,
                 "movl $0x3f800000, 12(%%esp)\n\t"// since the xmm4 registers allow
                 "movaps (%%esp), %%xmm3\n\t"     // us to operate on 4 floats at
                 "movaps (%%esp), %%xmm4\n\t"     // a time

                 "addl $16, %%esp\n\t"            // deallocate 4 floats on stack

                 "xorps %%xmm5, %%xmm5\n\t"       // set xmm5 to 0

                 // Loop over individual floats till we reach the right alignment
                 // that was computed in n0. If we don't start handling 4 floats
                 // at a time with SSE on a 4 bytes boundary, we get a crash
                 // in movaps (here, we use only movss, moving only 1 float at a
                 // time).
                 "0:\n\t"
                 "test %%ecx, %%ecx\n\t"          // if n0 == 0, jump to next loop
                 "jz 1f\n\t"

                 "movss (%%esi), %%xmm0\n\t"      // move a single float to xmm0
                 "cmpss $1, %%xmm0, %%xmm1\n\t"   // compare to threshold
                 "andps %%xmm1, %%xmm3\n\t"       // and with all 1s
                 "addss %%xmm3, %%xmm5\n\t"       // add result to xmm5 (=count!)
                 "movaps %%xmm2, %%xmm1\n\t"      // restore threshold in xmm1
                 "movaps %%xmm4, %%xmm3\n\t"      // restore all 1s in xmm3
                 "addl $4, %%esi\n\t"             // move to next float (4 bytes)
                 "decl %%ecx\n\t"                 // decrement ecx, which started at n0
                 "ja 0b\n\t"                      // jump if not done yet

                 // Loop over 4 floats at a time: this time, we have reached
                 // the proper alignment for movaps, so we can operate in parallel
                 // on 4 floats at a time. The code is the same as the previous loop
                 // except that the "ss" instructions are now "ps" instructions.
                 "1:\n\t"
                 "test %%edx, %%edx\n\t"
                 "jz 2f\n\t"

                 "movaps (%%esi), %%xmm0\n\t"     // note movaps, not movss
                 "cmpps $1, %%xmm0, %%xmm1\n\t"
                 "andps %%xmm1, %%xmm3\n\t"
                 "addps %%xmm3, %%xmm5\n\t"       // addps, not addss
                 "movaps %%xmm2, %%xmm1\n\t"
                 "movaps %%xmm4, %%xmm3\n\t"
                 "addl $16, %%esi\n\t"            // jump over 4 floats
                 "subl $4, %%edx\n\t"             // decrement edx (n1) by 4
                 "ja 1b\n\t"

                 // Tally up count so far into last float of xmm5: we were
                 // doing operations in parallels on the 4 floats in the xmm
                 // registers, resulting in 4 partial counts in xmm5.
                 "xorps %%xmm0, %%xmm0\n\t"
                 "haddps %%xmm0, %%xmm5\n\t"
                 "haddps %%xmm0, %%xmm5\n\t"

                 // Last loop, for stragglers in case the array is not evenly
                 // divisible by 4. We are back to operating on a single float
                 // at a time, using movss and addss.
                 "2:\n\t"
                 "test %%edi, %%edi\n\t"
                 "jz 3f\n\t"

                 "movss (%%esi), %%xmm0\n\t"
                 "cmpss $1, %%xmm0, %%xmm1\n\t"
                 "andps %%xmm1, %%xmm3\n\t"
                 "addss %%xmm3, %%xmm5\n\t"
                 "movaps %%xmm2, %%xmm1\n\t"
                 "movaps %%xmm4, %%xmm3\n\t"
                 "addl $4, %%esi\n\t"
                 "decl %%edi\n\t"
                 "ja 0b\n\t"

                 // Push result from xmm5 to variable count in memory.
                 "3:\n\t"
                 "movss %%xmm5, %0\n\t"

                 : "=m" (count)
                 : "S" (begin), "a" (threshold), "c" (n0), "d" (n1), "D" (n2)
                 :
                 );

    return (int) count;

  } else {
    return std::count_if(begin, end, std::bind2nd(std::greater<float>(), threshold));
  }
#else
  return std::count_if(begin, end, std::bind2nd(std::greater<float>(), threshold));
#endif
}

//--------------------------------------------------------------------------------
/**
 * Counts the number of non-zeros in a vector.
 * Doesn't work with vector<bool>.
 */
template<typename T>
inline size_t count_non_zeros(const std::vector<T>& x) {
  assert(sizeof(T) == 4); // because calls asm: T != size_t (8 bytes)
  float *begin = (float *) &x[0];
  float *end = begin + x.size();
  return count_gt(begin, end, 0);
}

//--------------------------------------------------------------------------------
template<typename T>
inline size_t n_zeros(const std::vector<T>& x) {
  return count_non_zeros(x);
}

//--------------------------------------------------------------------------------
/**
 * TODO: Use SSE. Maybe requires having our own vector<bool> so that we can avoid
 * the shenanigans with the bit references and iterators.
 */
template<>
inline size_t count_non_zeros(const std::vector<bool>& x) {
  size_t count = 0;
  for (size_t i = 0; i != x.size(); ++i)
    count += x[i];
  return count;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline size_t count_non_zeros(const std::vector<std::pair<T1, T2> >& x) {
  size_t count = 0;
  for (size_t i = 0; i != x.size(); ++i)
    if (!is_zero(x[i]))
      ++count;
  return count;
}

//--------------------------------------------------------------------------------
/**
 * Counts the number of values less than a given threshold in a given range.
 */
template<typename It>
inline size_t
count_lt(It begin, It end, const typename std::iterator_traits<It>::value_type& thres) {
  typedef typename std::iterator_traits<It>::value_type value_type;
  return std::count_if(begin, end, std::bind2nd(std::less<value_type>(), thres));
}

//--------------------------------------------------------------------------------
// Rounding
//--------------------------------------------------------------------------------
/**
 */
template<typename It>
inline void
round_01(It begin, It end,
         const typename std::iterator_traits<It>::value_type& threshold = .5) {
  {
    assert(begin <= end);
  }

  typename std::iterator_traits<It>::value_type val;

  while (begin != end) {
    val = *begin;
    if (val >= threshold)
      val = 1;
    else
      val = 0;
    *begin = val;
    ++begin;
  }
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T>
inline void round_01(T& a, const typename T::value_type& threshold = .5) {
  round_01(a.begin(), a.end(), threshold);
}

//--------------------------------------------------------------------------------
// Addition...
//--------------------------------------------------------------------------------
/**
 * Computes the sum of the elements in a range.
 * vDSP is much faster than C++, even optimized by gcc, but for now this works
 * only with float (rather than double), and only on darwin86. With these
 * restrictions the speed-up is usually better than 5X over optimized C++.
 * vDSP also handles unaligned vectors correctly, and has good performance
 * also when the vectors are small, not just when they are big.
 */
inline float sum(float *begin, float *end) {
  {
    assert(begin <= end);
  }

#ifdef PLATFORM_darwin86

  float result = 0;
  vDSP_sve(begin, 1, &result, (end - begin));
  return result;

#else

  float result = 0;
  for (; begin != end; ++begin)
    result += *begin;
  return result;

#endif
}

//--------------------------------------------------------------------------------
/**
 * Compute the sum of a whole container.
 * Here we revert to C++, which is going to be slower than the preceding function,
 * but it will work for a container of anything, that container not necessarily
 * being a contiguous vector of numbers.
 */
template<typename T>
inline typename T::value_type sum(const T& x) {
  typename T::value_type result = 0;
  typename T::const_iterator it;
  for (it = x.begin(); it != x.end(); ++it)
    result += *it;
  return result;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2, typename T3>
inline void sum(const std::vector<T1>& a, const std::vector<T2>& b,
                size_t begin, size_t end, std::vector<T3>& c) {
  for (size_t i = begin; i != end; ++i)
    c[i] = a[i] + b[i];
}

//--------------------------------------------------------------------------------
/**
 * Computes the product of the elements in a range.
 */
template<typename It>
inline typename std::iterator_traits<It>::value_type product(It begin, It end) {
  {
    assert(begin <= end);
  }

  typename std::iterator_traits<It>::value_type p(1);

  for (; begin != end; ++begin)
    p *= *begin;

  return p;
}

//--------------------------------------------------------------------------------
/**
 * Computes the product of all the elements in a container.
 */
template<typename T>
inline typename T::value_type product(const T& x) {
  return product(x.begin(), x.end());
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It>
inline void add_val(It begin, It end,
                    const typename std::iterator_traits<It>::value_type& val) {
  {
    assert(begin <= end);
  }

  for (; begin != end; ++begin)
    *begin += val;
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T>
inline void add_val(T& x, const typename T::value_type& val) {
  add_val(x.begin(), x.end(), val);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It>
inline void subtract_val(It begin, It end,
                         const typename std::iterator_traits<It>::value_type& val) {
  {
    assert(begin <= end);
  }

  for (; begin != end; ++begin)
    *begin -= val;
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T>
inline void subtract_val(T& x, const typename T::value_type& val) {
  subtract_val(x.begin(), x.end(), val);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It>
inline void negate(It begin, It end) {
  {
    assert(begin <= end);
  }

  for (; begin != end; ++begin)
    *begin = -*begin;
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T>
inline void negate(T& x) {
  typename T::iterator it;
  for (it = x.begin(); it != x.end(); ++it)
    *it = -*it;
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It>
inline void
multiply_val(It begin, It end,
             const typename std::iterator_traits<It>::value_type& val) {
  {
    assert(begin <= end);
  }

  for (; begin != end; ++begin)
    *begin *= val;
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T>
inline void multiply_val(T& x, const typename T::value_type& val) {
  multiply_val(x.begin(), x.end(), val);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It>
inline void
divide_val(It begin, It end,
           const typename std::iterator_traits<It>::value_type& val) {
  {
    assert(begin <= end);
    assert(val != 0);
  }

  for (; begin != end; ++begin)
    *begin /= val;
}

//--------------------------------------------------------------------------------
// TODO: what if val == 0?
/**
 */
template<typename T>
inline void divide_val(T& x, const typename T::value_type& val) {
  divide_val(x.begin(), x.end(), val);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It1, typename It2>
inline void add(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    assert(begin1 <= end1);
    assert(end1 - begin1 <= end2 - begin2);
  }

  for (; begin1 != end1; ++begin1, ++begin2)
    *begin1 += *begin2;
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T1, typename T2>
inline void add(T1& x, const T2& y) {
  add(x.begin(), x.end(), y.begin(), y.end());
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It1, typename It2>
inline void subtract(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    assert(begin1 <= end1);
    assert(end1 - begin1 <= end2 - begin2);
  }

  for (; begin1 != end1; ++begin1, ++begin2)
    *begin1 -= *begin2;
}

//--------------------------------------------------------------------------------
// TODO: should we have the same argument ordering as copy??
/**
 */
template<typename T1, typename T2>
inline void subtract(T1& x, const T2& y) {
  subtract(x.begin(), x.end(), y.begin(), y.end());
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It1, typename It2>
inline void multiply(It1 begin1, It1 end1, It2 begin2, It2 end2) {
  {
    assert(begin1 <= end1);
    assert(end1 - begin1 <= end2 - begin2);
  }

  for (; begin1 != end1; ++begin1, ++begin2)
    *begin1 *= *begin2;
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T1, typename T2>
inline void multiply(T1& x, const T2& y) {
  multiply(x.begin(), x.end(), y.begin(), y.end());
}

//--------------------------------------------------------------------------------
template<typename It1, typename It2, typename It3>
inline void multiply(It1 begin1, It1 end1, It2 begin2, It2 end2,
                     It3 begin3, It3 end3) {
  {
    assert(begin1 <= end1);
    assert(end1 - begin1 <= end2 - begin2);
    assert(end1 - begin1 <= end3 - begin3);
  }

  typedef typename std::iterator_traits<It3>::value_type value_type;

  for (; begin1 != end1; ++begin1, ++begin2, ++begin3)
    *begin3 = (value_type) *begin1 * *begin2;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2, typename T3>
inline void multiply(const T1& x, const T2& y, T3& z) {
  multiply(x.begin(), x.end(), y.begin(), y.end(), z.begin(), z.end());
}

//--------------------------------------------------------------------------------
/**
 * Given a vector of pairs <index, value> and a value val, multiplies the values
 * by val, but only if index is in indices. Needs x and indices to be sorted
 * in order of increasing indices.
 */
template<typename I, typename T>
inline void
multiply_val(T val, const std::vector<I>& indices, SparseVector<I, T>& x) {
  I n1 = indices.size(), n2 = x.size(), i1 = 0, i2 = 0;

  while (i1 != n1 && i2 != n2)
    if (x[i2].first < indices[i1]) {
      ++i2;
    } else if (indices[i1] < x[i2].first) {
      ++i1;
    } else {
      x[i2].second *= val;
      ++i1;
      ++i2;
    }
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It1, typename It2>
inline void divide(It1 begin1, It1 end1, It2 begin2, It2 end2,
                   typename std::iterator_traits<It1>::value_type fuzz = 0) {
  {
    assert(begin1 <= end1);
    assert(end1 - begin1 <= end2 - begin2);
  }

  if (fuzz == 0)
    for (; begin1 != end1; ++begin1, ++begin2)
      *begin1 /= *begin2;
  else
    for (; begin1 != end1; ++begin1, ++begin2)
      *begin1 /= (*begin2 + fuzz);
}

//--------------------------------------------------------------------------------
// What if y contains one or more zeros?
/**
 */
template<typename T1, typename T2>
inline void divide(T1& x, const T2& y, typename T1::value_type fuzz = 0) {
  divide(x.begin(), x.end(), y.begin(), y.end(), fuzz);
}

//--------------------------------------------------------------------------------
template<typename It1>
inline void divide_by_max(It1 begin, It1 end) {
  {
    assert(begin <= end);
  }

  typename std::iterator_traits<It1>::value_type max_val =
  *(std::max_element(begin, end));

  if (!nearlyZero(max_val))
    for (It1 it = begin; it != end; ++it)
      *it /= max_val;
}

//--------------------------------------------------------------------------------
template<typename T1>
inline void divide_by_max(T1& v) {
  divide_by_max(v.begin(), v.end());
}

//--------------------------------------------------------------------------------
/**
 */
template<typename It1, typename It2, typename TFuncIsNearlyZero, typename TFuncHandleZero>
inline void inverseNZ(It1 begin1, It1 end1, It2 out, It2 out_end,
                      TFuncIsNearlyZero fIsZero, TFuncHandleZero fHandleZero) {
  {
    assert(begin1 <= end1);
    assert(out <= out_end);
    assert(end1 - begin1 == out_end - out);
  }

  const typename std::iterator_traits<It2>::value_type one(1.0);

  for (; begin1 != end1; ++begin1, ++out) {
    if (fIsZero(*begin1))
      *out = fHandleZero(*begin1); // Can't pass one?
    else
      *out = one / *begin1;
  }
}

//--------------------------------------------------------------------------------
/**
 * Computes the reciprocal of each element of vector 'x' and writes
 * result into vector 'out'.
 * 'out' must be of at least the size of 'x'.
 * Does not resize 'out'; behavior is undefined if 'out' is not of
 * the correct size.
 * Uses only 'value_type', 'begin()' and 'end()' of 'x' and
 * 'value_type' and 'begin()' of 'out'.
 * Checks the value of each element of 'x' with 'fIsNearlyZero', and
 * if 'fIsNearlyZero' returns false, computes the reciprocal as
 * T2::value_type(1.0) / element value.
 * If 'fIsNearlyZero' returns true, computes uses the output of
 * 'fHandleZero(element value)' as the result for that element.
 *
 * Usage: inverseNZ(input, output,
 *            IsNearlyZero< DistanceToZero<double> >(),
 *            Identity<double>());
 */
template<typename T1, typename T2, typename TFuncIsNearlyZero, typename TFuncHandleZero>
inline void inverseNZ(const T1& x, T2& out,
                      TFuncIsNearlyZero fIsNearlyZero, TFuncHandleZero fHandleZero) {
  inverseNZ(x.begin(), x.end(), out.begin(), out.end(), fIsNearlyZero, fHandleZero);
}

//--------------------------------------------------------------------------------
template<typename It1, typename It2>
inline void inverse(It1 begin1, It1 end1, It2 out, It2 out_end,
                    const typename std::iterator_traits<It2>::value_type one = 1.0) {
  {
    assert(begin1 <= end1);
    assert(out <= out_end);
    assert(end1 - begin1 == out_end - out);
  }

  for (; begin1 != end1; ++begin1, ++out)
    *out = one / *begin1;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline void inverse(const T1& x, T2& out, const typename T2::value_type one = 1.0) {
  inverse(x.begin(), x.end(), out.begin(), out.end(), one);
}

//--------------------------------------------------------------------------------
/**
 * x += k y
 */
template<typename It1, typename It2>
inline void add_ky(const typename std::iterator_traits<It1>::value_type& k,
                   It1 y, It1 y_end, It2 x, It2 x_end) {
  {
    assert(y <= y_end);
    assert(x <= x_end);
    assert(y_end - y <= x - x_end);
  }

  while (y != y_end) {
    *x += k * *y;
    ++x;
    ++y;
  }
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T1, typename T2>
inline void add_ky(const typename T1::value_type& k, const T2& y, T1& x) {
  add_ky(k, y.begin(), y.end(), x.begin(), x.end());
}

//--------------------------------------------------------------------------------
/**
 * x2 = x1 + k y
 */
template<typename It1, typename It2, typename It3>
inline void add_ky(It1 x1, It1 x1_end,
                   const typename std::iterator_traits<It1>::value_type& k,
                   It2 y, It3 x2) {
  while (x1 != x1_end) {
    *x2 = *x1 + k * *y;
    ++x2;
    ++x1;
    ++y;
  }
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T1, typename T2, typename T3>
inline void add_ky(const T1& x1,
                   const typename T1::value_type& k, const T2& y,
                   T3& x2) {
  ////assert(y.size() >= x.size());

  add_ky(x1.begin(), x1.end(), k, y.begin(), x2.begin());
}

// TODO: write binary operations x = y + z ...

//--------------------------------------------------------------------------------
/**
 * x = a * x + y
 *
 * TODO: write the rest of BLAS level 1
 */
template<typename T1, typename T2>
inline void axpy(T1& x, const typename T1::value_type& a, const T2& y) {
  ////assert(y.size() >= x.size());

  typename T1::iterator it_x = x.begin(), it_x_end = x.end();
  typename T2::const_iterator it_y = y.begin();

  while (it_x != it_x_end) {
    *it_x = a * *it_x + *it_y;
    ++it_x;
    ++it_y;
  }
}

//--------------------------------------------------------------------------------
/**
 * x = a * x + b * y
 */
template<typename X, typename Y>
inline void axby(const typename std::iterator_traits<X>::value_type& a,
                 X x, X x_end,
                 const typename std::iterator_traits<X>::value_type& b,
                 Y y) {
  while (x != x_end) {
    *x = a * *x + b * *y;
    ++x;
    ++y;
  }
}

//--------------------------------------------------------------------------------
/**
 */
template<typename T1, typename T2>
inline void axby(const typename T1::value_type& a, T1& x,
                 const typename T1::value_type& b, const T2& y) {
  ////assert(y.size() >= x.size());

  axby(a, x.begin(), x.end(), b, y.begin());
}

//--------------------------------------------------------------------------------
/**
 * exp(k * x) for all the elements of a range.
 */
template<typename It>
inline void range_exp(typename std::iterator_traits<It>::value_type k,
                      It begin, It end) {
  typedef typename std::iterator_traits<It>::value_type value_type;

  Exp <value_type> e_f;

  for (; begin != end; ++begin)
    *begin = e_f(k * *begin);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename C>
inline void range_exp(typename C::value_type k, C& c) {
  range_exp(k, c.begin(), c.end());
}

//--------------------------------------------------------------------------------
/**
 * k1 * exp(k2 * x) for all the elements of a range.
 */
template<typename It>
inline void range_exp(typename std::iterator_traits<It>::value_type k1,
                      typename std::iterator_traits<It>::value_type k2,
                      It begin, It end) {
  typedef typename std::iterator_traits<It>::value_type value_type;

  Exp <value_type> e_f;

  for (; begin != end; ++begin)
    *begin = k1 * e_f(k2 * *begin);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename C>
inline void range_exp(typename C::value_type k1, typename C::value_type k2, C& c) {
  range_exp(k1, k2, c.begin(), c.end());
}

//--------------------------------------------------------------------------------
// Inner product
//--------------------------------------------------------------------------------
/**
 * Bypasses the STL API and its init value.
 * TODO: when range is empty??
 */
template<typename It1, typename It2>
inline typename std::iterator_traits<It1>::value_type
inner_product(It1 it_x, It1 it_x_end, It2 it_y) {
  typename std::iterator_traits<It1>::value_type n(0);

  while (it_x != it_x_end) {
    n += *it_x * *it_y;
    ++it_x;
    ++it_y;
  }

  return n;
}

//--------------------------------------------------------------------------------
/**
 * In place transform of a range.
 */
template<typename F1, typename It>
inline void transform(It begin, It end, F1 f) {
  {
    assert(begin <= end);
  }

  for (; begin != end; ++begin)
    *begin = f(*begin);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename F1, typename T1>
inline void transform(T1& a, F1 f) {
  typename T1::iterator ia = a.begin(), iae = a.end();

  for (; ia != iae; ++ia)
    *ia = f(*ia);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename F1, typename T1, typename T2>
inline void transform(const T1& a, T2& b, F1 f) {
  ////assert(b.size() >= a.size());

  typename T1::const_iterator ia = a.begin(), iae = a.end();
  typename T2::iterator ib = b.begin();

  for (; ia != iae; ++ia, ++ib)
    *ib = f(*ia);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename F2, typename T1, typename T2, typename T3>
inline void transform(const T1& a, const T2& b, T3& c, F2 f) {
  ////assert(c.size() >= a.size());
  ////assert(b.size() >= a.size());
  ////assert(c.size() >= b.size());

  typename T1::const_iterator ia = a.begin(), iae = a.end();
  typename T2::const_iterator ib = b.begin();
  typename T3::iterator ic = c.begin();

  for (; ia != iae; ++ia, ++ib, ++ic)
    *ic = f(*ia, *ib);
}

//--------------------------------------------------------------------------------
/**
 */
template<typename F3, typename T1, typename T2, typename T3, typename T4>
inline void transform(const T1& a, const T2& b, const T3& c, T4& d, F3 f) {
  ////assert(d.size() >= a.size());
  ////assert(d.size() >= b.size());
  ////assert(d.size() >= c.size());
  ////assert(b.size() >= a.size());
  ////assert(c.size() >= a.size());

  typename T1::const_iterator ia = a.begin(), iae = a.end();
  typename T2::const_iterator ib = b.begin();
  typename T3::const_iterator ic = c.begin();
  typename T4::iterator id = d.begin();

  for (; ia != iae; ++ia, ++ib, ++ic, ++id)
    *id = f(*ia, *ib, *ic);
}

//--------------------------------------------------------------------------------
/*
 * Removes consecutive duplicates in a vector. Requires equality operator to be
 * defined on T. Keeps the _first_ element in a series of duplicates.
 */
template<typename T>
inline void remove_consecutive_duplicates(std::vector<T>& x) {
  size_t i1 = 0, i2 = 0;
  while (i1 < x.size()) {
    x[i2++] = x[i1++];
    while (i1 < x.size() && x[i1] == x[i2 - 1])
      ++i1;
  }

  x.resize(i2);
}

//--------------------------------------------------------------------------------
/*
 * Removes consecutive duplicates in a list. Requires equality operator to be
 * defined on T. Keeps the _first_ element in a series of duplicates.
 */
template<typename T>
inline void remove_consecutive_duplicates(std::list<T>& x) {
  typename std::list<T>::iterator prev, it;

  it = x.begin();
  while (it != x.end()) {
    prev = it;
    ++it;
    while (it != x.end() && *prev == *it)
      it = x.erase(it);
  }
}

//--------------------------------------------------------------------------------
// min_element / max_element
//--------------------------------------------------------------------------------
template<typename T>
inline T max_value(const std::vector<T>& x) {
  return *std::max_element(x.begin(), x.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline T min_value(const std::vector<T>& x) {
  return *std::min_element(x.begin(), x.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline T range(const std::vector<T>& x) {
  T min_v = std::numeric_limits<T>::max();
  T max_v = -std::numeric_limits<T>::max();

  T *p = &x[0], *p_end = p + x.size();
  for (; p != p_end; ++p) {
    if (*p < min_v)
      min_v = *p;
    if (*p > max_v)
      max_v = *p;
  }

  return max_v - min_v;
}

//--------------------------------------------------------------------------------
/**
 * Returns the position at which f takes its minimum between first and last.
 */
template<typename ForwardIterator, typename F>
inline ForwardIterator
min_element(ForwardIterator first, ForwardIterator last, F f) {
  {
    assert(first <= last);
  }

  typedef typename ForwardIterator::value_type value_type;

  ForwardIterator min_it = first;
  value_type min_val = f(*first);

  while (first != last) {
    value_type val = f(*first);
    if (val < min_val) {
      min_it = first;
      min_val = val;
    }
    ++first;
  }

  return min_it;
}

//--------------------------------------------------------------------------------
/**
 * Returns the position at which f takes its maximum between first and last.
 */
template<typename ForwardIterator, typename F>
inline ForwardIterator
max_element(ForwardIterator first, ForwardIterator last, F f) {
  {
    assert(first <= last);
  }

  typedef typename ForwardIterator::value_type value_type;

  ForwardIterator max_it = first;
  value_type max_val = f(*first);

  while (first != last) {
    value_type val = f(*first);
    if (val > max_val) {
      max_it = first;
      max_val = val;
    }
    ++first;
  }

  return max_it;
}

//--------------------------------------------------------------------------------
/**
 * Finds the min element in a container.
 */
template<typename C>
inline size_t arg_min(const C& c) {
  if (c.empty())
    return (size_t) 0;
  else
    return (size_t) (std::min_element(c.begin(), c.end()) - c.begin());
}

//--------------------------------------------------------------------------------
/**
 * Finds the maximum element in a container.
 */
template<typename C>
inline size_t arg_max(const C& c) {
  if (c.empty())
    return (size_t) 0;
  else
    return (size_t) (std::max_element(c.begin(), c.end()) - c.begin());
}

//--------------------------------------------------------------------------------
/**
 * Writes the component-wise minimum to the output vector.
 */
template<typename It1, typename It2, typename It3>
inline void minimum(It1 begin1, It1 end1, It2 begin2, It3 out) {
  {
    assert(begin1 <= end1);
  }

  typedef typename std::iterator_traits<It3>::value_type T;
  for (; begin1 != end1; ++begin1, ++begin2, ++out) {
    *out = std::min<T>(*begin1, *begin2);
  }
}

//--------------------------------------------------------------------------------
/**
 * Writes the component-wise minimum to the output vector.
 */
template<typename T1, typename T2, typename T3>
inline void minimum(const T1& x, const T2& y, T3& out) {
  minimum(x.begin(), x.end(), y.begin(), out.begin());
}

//--------------------------------------------------------------------------------
// contains
//--------------------------------------------------------------------------------
template<typename C>
inline bool contains(const C& c, typename C::value_type& v) {
  return std::find(c.begin(), c.end(), v) != c.end();
}

//--------------------------------------------------------------------------------
template<typename C1, typename C2>
inline bool is_subsequence(const C1& seq, const C2& sub) {
  return std::search(seq.begin(), seq.end(), sub.begin(), sub.end()) != seq.end();
}

//--------------------------------------------------------------------------------
template<typename C1, typename C2>
inline bool is_subsequence_of(const C1& c, const C2& sub) {
  bool found = false;
  typename C1::const_iterator it, end = c.end();
  for (it = c.begin(); it != end; ++it)
    if (is_subsequence(*it, sub))
      found = true;
  return found;
}

//--------------------------------------------------------------------------------
// sample
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
// DENSE LOGICAL AND/OR
//--------------------------------------------------------------------------------
/**
 * For each corresponding elements of x and y, put the logical and of those two
 * elements at the corresponding position in z. This is faster than the numpy
 * logical_and, which doesn't seem to be using SSE.
 *
 * x, y and z are arrays of floats, but with 0/1 values.
 *
 * If any of the vectors is not aligned on a 16 bytes boundary, the function
 * reverts to slow C++. This can happen when using it with slices of numpy
 * arrays.
 *
 * Doesn't work on 64 bits platforms, doesn't work on win32.
 *
 * TODO: find 16 bytes aligned block that can be sent to SSE.
 * TODO: support other platforms than just darwin86 for the fast path.
 */
template<typename InputIterator, typename OutputIterator>
inline void logical_and(InputIterator x, InputIterator x_end,
                        InputIterator y, InputIterator y_end,
                        OutputIterator z, OutputIterator z_end) {
  {
    assert(x_end - x == y_end - y);
    assert(x_end - x == z_end - z);
  }

  // See comments in count_gt. We need both conditional compilation and
  // SSE_LEVEL check.
#ifdef PLATFORM_darwin86

  if (SSE_LEVEL >= 3) {

    // n is the total number of floats to process
    // n1 is 0 if any of the arrays x,y,z is not aligned on a 4 bytes
    // boundary, or the number of floats we'll be able to process in parallel
    // using the xmm.
    int n = (int)(x_end - x);
    int n1 = 0;
    if (((long)x) % 16 == 0 && ((long)y) % 16 == 0 && ((long)z) % 16 == 0)
      n1 = 16 * (n / 16);

    // If we are not aligned on 4 bytes, n1 == 0, and we simply
    // skip the asm.
    if (n1 > 0) {

      asm volatile(
                   "pusha\n\t"                   // save all registers

                   "0:\n\t"
                   "movaps (%%esi), %%xmm0\n\t"  // move 4 floats of x to xmm0
                   "andps (%%edi), %%xmm0\n\t"   // parallel and with 4 floats of y
                   "movaps 16(%%esi), %%xmm1\n\t"// play again with next 4 floats
                   "andps 16(%%edi), %%xmm1\n\t"
                   "movaps 32(%%esi), %%xmm2\n\t"// and next 4 floats
                   "andps 32(%%edi), %%xmm2\n\t"
                   "movaps 48(%%esi), %%xmm3\n\t"// and next 4 floats: we've and'ed
                   "andps 48(%%edi), %%xmm3\n\t" // 16 floats of x and y at this point

                   "movaps %%xmm0, (%%ecx)\n\t"  // simply move 4 floats at a time to z
                   "movaps %%xmm1, 16(%%ecx)\n\t"// and next 4 floats
                   "movaps %%xmm2, 32(%%ecx)\n\t"// and next 4 floats
                   "movaps %%xmm3, 48(%%ecx)\n\t"// and next 4: moved 16 floats to z

                   "addl $64, %%esi\n\t"         // increment pointer into x by 16 floats
                   "addl $64, %%edi\n\t"         // increment pointer into y
                   "addl $64, %%ecx\n\t"         // increment pointer into z
                   "subl $16, %%edx\n\t"         // we've processed 16 floats
                   "ja 0b\n\t"                   // loop

                   "popa\n\t"                    // restore registers

                   :
                   : "S" (x), "D" (y), "c" (z), "d" (n1)
                   :
                   );
    }

    // Finish up for stragglers in case the array length was not
    // evenly divisible by 4
    for (int i = n1; i != n; ++i)
      *(z+i) = *(x+i) && *(y+i);

  } else {

    for (; x != x_end; ++x, ++y, ++z)
      *z = (*x) && (*y);

  }
#else
  for (; x != x_end; ++x, ++y, ++z)
    *z = (*x) && (*y);
#endif
}

//--------------------------------------------------------------------------------
/**
 * Same as previous logical_and, but puts the result back into y.
 * Same comments.
 */
template<typename Iterator>
inline void in_place_logical_and(Iterator x, Iterator x_end,
                                 Iterator y, Iterator y_end) {
  {
    assert(x_end - x == y_end - y);
  }

  // See comments in count_gt. We need conditional compilation
  // _AND_ SSE_LEVEL check.
#ifdef PLATFORM_darwin86

  if (SSE_LEVEL >= 3) {

    // See comments in logical_and.
    int n = (int)(x_end - x);
    int n1 = 0;
    if (((long)x) % 16 == 0 && ((long)y) % 16 == 0)
      n1 = 16 * (n / 16);

    if (n1 > 0) {

      asm volatile(
                   "pusha\n\t"

                   "0:\n\t"
                   "movaps (%%esi), %%xmm0\n\t"
                   "movaps 16(%%esi), %%xmm1\n\t"
                   "movaps 32(%%esi), %%xmm2\n\t"
                   "movaps 48(%%esi), %%xmm3\n\t"
                   "andps (%%edi), %%xmm0\n\t"
                   "andps 16(%%edi), %%xmm1\n\t"
                   "andps 32(%%edi), %%xmm2\n\t"
                   "andps 48(%%edi), %%xmm3\n\t"
                   "movaps %%xmm0, (%%edi)\n\t"
                   "movaps %%xmm1, 16(%%edi)\n\t"
                   "movaps %%xmm2, 32(%%edi)\n\t"
                   "movaps %%xmm3, 48(%%edi)\n\t"

                   "addl $64, %%esi\n\t"
                   "addl $64, %%edi\n\t"
                   "subl $16, %%edx\n\t"
                   "prefetch (%%esi)\n\t"
                   "ja 0b\n\t"

                   "popa\n\t"

                   :
                   : "S" (x), "D" (y), "d" (n1)
                   :
                   );
    }

    for (int i = n1; i != n; ++i)
      *(y+i) = *(x+i) && *(y+i);

  } else {

    for (; x != x_end; ++x, ++y)
      *y = (*x) && *(y);
  }
#else
  for (; x != x_end; ++x, ++y)
    *y = (*x) && *(y);
#endif
}

//--------------------------------------------------------------------------------
/**
 * A specialization tuned for unsigned char.
 * TODO: keep only one code that computes the right offsets based on
 * the iterator value type?
 * TODO: vectorize, but watch out for alignments
 */
inline void in_place_logical_and(const ByteVector& x, ByteVector& y,
                                 int begin = -1, int end = -1) {
  if (begin == -1)
    begin = 0;

  if (end == -1)
    end = (int) x.size();

  for (int i = begin; i != end; ++i)
    y[i] &= x[i];
}

//--------------------------------------------------------------------------------
/**
 * TODO: write with SSE for big enough vectors.
 */
inline void in_place_logical_or(const ByteVector& x, ByteVector& y,
                                int begin = -1, int end = -1) {
  if (begin == -1)
    begin = 0;

  if (end == -1)
    end = (int) x.size();

  for (int i = begin; i != end; ++i)
    y[i] |= x[i];
}

//--------------------------------------------------------------------------------
inline void
logical_or(size_t n, const ByteVector& x, const ByteVector& y, ByteVector& z) {
  for (size_t i = 0; i != n; ++i)
    z[i] = x[i] || y[i];
}

//--------------------------------------------------------------------------------
inline void in_place_logical_or(size_t n, const ByteVector& x, ByteVector& y) {
  for (size_t i = 0; i != n; ++i)
    y[i] |= x[i];
}

//--------------------------------------------------------------------------------
// SPARSE OR/AND
//--------------------------------------------------------------------------------
template<typename InputIterator1, typename InputIterator2, typename OutputIterator>
inline size_t sparse_or(InputIterator1 begin1, InputIterator1 end1,
                        InputIterator2 begin2, InputIterator2 end2,
                        OutputIterator out, OutputIterator out_end) {
  { // Pre-conditions
    assert(0 <= end1 - begin1);
    assert(0 <= end2 - begin2);
    assert(0 <= out_end - out);
  } // End pre-conditions

  typedef typename std::iterator_traits<OutputIterator>::value_type value_type;

  OutputIterator out_begin = out;

  while (begin1 != end1 && begin2 != end2) {

    if (*begin1 < *begin2) {
      *out++ = (value_type) *begin1++;
    } else if (*begin2 < *begin1) {
      *out++ = (value_type) *begin2++;
    } else {
      *out++ = (value_type) *begin1++;
      ++begin2;
    }
  }

  for (; begin1 != end1; ++begin1)
    *out++ = (value_type) *begin1;

  for (; begin2 != end2; ++begin2)
    *out++ = (value_type) *begin2;

  return (size_t) (out - out_begin);
}

//--------------------------------------------------------------------------------
template<typename U>
inline size_t sparse_or(const std::vector<U>& x1, const std::vector<U>& x2,
                        std::vector<U>& out) {
  return sparse_or(x1.begin(), x1.end(), x2.begin(), x2.end(),
                   out.begin(), out.end());
}

//--------------------------------------------------------------------------------
template<typename InputIterator1, typename InputIterator2, typename OutputIterator>
inline size_t sparse_and(InputIterator1 begin1, InputIterator1 end1,
                         InputIterator2 begin2, InputIterator2 end2,
                         OutputIterator out, OutputIterator out_end) {
  { // Pre-conditions
    assert(0 <= end1 - begin1);
    assert(0 <= end2 - begin2);
    assert(0 <= out_end - out);
  } // End pre-conditions

  typedef typename std::iterator_traits<OutputIterator>::value_type value_type;

  OutputIterator out_begin = out;

  while (begin1 != end1 && begin2 != end2) {

    if (*begin1 < *begin2) {
      ++begin1;
    } else if (*begin2 < *begin1) {
      ++begin2;
    } else {
      *out++ = (value_type) *begin1++;
      ++begin2;
    }
  }

  return (size_t) (out - out_begin);
}

//--------------------------------------------------------------------------------
template<typename U>
inline size_t sparse_and(const std::vector<U>& x1, const std::vector<U>& x2,
                         std::vector<U>& out) {
  return sparse_and(x1.begin(), x1.end(), x2.begin(), x2.end(),
                    out.begin(), out.end());
}

//--------------------------------------------------------------------------------
// sort, partial_sort
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
/*
template <typename C>
inline void sort(C& c)
{
  std::sort(c.begin(), c.end());
}
*/

//--------------------------------------------------------------------------------
template<typename C, typename F>
inline void sort(C& c, F f) {
  std::sort(c.begin(), c.end(), f);
}

//--------------------------------------------------------------------------------
template<typename It>
inline void sort_on_first(It x_begin, It x_end, int direction = 1) {
  typedef typename std::iterator_traits<It>::value_type P;
  typedef typename P::first_type I;
  typedef typename P::second_type F;

  typedef select1st <std::pair<I, F>> sel1st;

  if (direction == -1) {
    std::sort(x_begin, x_end, predicate_compose<std::greater<I>, sel1st>());
  } else {
    std::sort(x_begin, x_end, predicate_compose<std::less<I>, sel1st>());
  }
}

//--------------------------------------------------------------------------------
template<typename It>
inline void sort_on_second(It x_begin, It x_end, int direction = 1) {
  typedef typename std::iterator_traits<It>::value_type P;
  typedef typename P::first_type I;
  typedef typename P::second_type F;

  typedef select2nd <std::pair<I, F>> sel2nd;

  if (direction == -1) {
    std::sort(x_begin, x_end, predicate_compose<std::greater<I>, sel2nd>());
  } else {
    std::sort(x_begin, x_end, predicate_compose<std::less<I>, sel2nd>());
  }
}

//--------------------------------------------------------------------------------
template<typename I, typename F>
inline void sort_on_first(size_t n, std::vector<std::pair<I, F> >& x,
                          int direction = 1) {
  sort_on_first(x.begin(), x.begin() + n, direction);
}

//--------------------------------------------------------------------------------
template<typename I, typename F>
inline void sort_on_first(std::vector<std::pair<I, F> >& x, int direction = 1) {
  sort_on_first(x.begin(), x.end(), direction);
}

//--------------------------------------------------------------------------------
template<typename I, typename F>
inline void sort_on_first(SparseVector<I, F>& x, int direction = 1) {
  sort_on_first(x.begin(), x.begin() + x.size(), direction);
}

//--------------------------------------------------------------------------------
template<typename I, typename C>
inline void partial_sort(I k, C& elts) {
  k = std::min(k, (I) elts.size());
  std::partial_sort(elts.begin(), elts.begin() + k, elts.end());
}

//--------------------------------------------------------------------------------
/**
 * Partial sort of a container given a functor.
 */
template<typename I, typename C, typename F>
inline void partial_sort(I k, C& elts, F f) {
  std::partial_sort(elts.begin(), elts.begin() + k, elts.end(), f);
}

//--------------------------------------------------------------------------------
/**
 * Partial sort of a range, that returns the values and the indices.
 */
template<typename size_type,
typename InputIterator,
typename OutputIterator,
typename Order>
inline void
partial_sort_2nd(size_type k,
                 InputIterator in_begin, InputIterator in_end,
                 OutputIterator out_begin, Order) {
  typedef typename std::iterator_traits<InputIterator>::value_type value_type;
  typedef select2nd <std::pair<size_type, value_type>> sel2nd;

  std::vector<std::pair<size_type, value_type> > v(in_end - in_begin);

  for (size_type i = 0; in_begin != in_end; ++in_begin, ++i)
    v[i] = std::make_pair(i, *in_begin);

  std::partial_sort(v.begin(), v.begin() + k, v.end(),
                    predicate_compose<Order, sel2nd>());

  for (size_type i = 0; i != k; ++i, ++out_begin)
    *out_begin = v[i];
}

//--------------------------------------------------------------------------------
/**
 * Partial sort of a container.
 */
template<typename C1, typename OutputIterator, typename Order>
inline void
partial_sort_2nd(size_t k, const C1& c1, OutputIterator out_begin, Order order) {
  partial_sort_2nd(k, c1.begin(), c1.end(), out_begin, order);
}

//--------------------------------------------------------------------------------
/**
 * Partial sort of a range of vectors, based on a given order predicate for
 * the vectors, putting the result into two iterators,
 * one for the indices and one for the element values.
 * Order needs to work for pairs (i.e., is a binary predicate).
 * start_offset specifies an optional for the indices that will be generated
 * for the pairs. This is useful when calling partial_sort_2nd repetitively
 * for different ranges inside a larger range.
 * If resort_on_first is true, the indices of the pairs are resorted,
 * otherwise, the indices might come out in any order.
 */
template<typename size_type,
typename InIter,
typename OutputIterator1,
typename OutputIterator2,
typename Order>
inline void
partial_sort(size_type k, InIter in_begin, InIter in_end,
             OutputIterator1 ind, OutputIterator2 nz,
             Order order, size_type start_offset = 0,
             bool resort_on_first = false) {
  typedef typename std::iterator_traits<InIter>::value_type value_type;
  typedef select1st <std::pair<size_type, value_type>> sel1st;

  std::vector<std::pair<size_type, value_type> > v(in_end - in_begin);

  for (size_type i = start_offset; in_begin != in_end; ++in_begin, ++i)
    v[i - start_offset] = std::make_pair(i, *in_begin);

  std::partial_sort(v.begin(), v.begin() + k, v.end(), order);

  if (resort_on_first) {
    std::sort(v.begin(), v.begin() + k,
              predicate_compose<std::less<size_type>, sel1st>());
  }

  for (size_type i = 0; i != k; ++i, ++ind, ++nz) {
    *ind = v[i].first;
    *nz = v[i].second;
  }
}

//--------------------------------------------------------------------------------
/**
 * A greater 2nd order that breaks ties, useful for debugging.
 */
template<typename T1, typename T2>
struct greater_2nd_no_ties
: public std::binary_function<bool, std::pair<T1, T2>, std::pair<T1, T2> > {
  inline bool operator()(const std::pair<T1, T2>& a, const std::pair<T1, T2>& b) const {
    if (a.second > b.second)
      return true;
    else if (a.second == b.second) if (a.first < b.first)
      return true;
    return false;
  }
};

//--------------------------------------------------------------------------------
/**
 * In place.
 */
template<typename I0, typename I, typename T>
inline void
partial_argsort(I0 k, SparseVector<I, T>& x, int direction = -1) {
  {
    assert(0 < k);
    assert(k <= x.size());
    assert(direction == -1 || direction == 1);
  }

  typedef I size_type;
  typedef T value_type;

  if (direction == -1) {

    greater_2nd_no_ties <size_type, value_type> order;
    std::partial_sort(x.begin(), x.begin() + k, x.begin() + x.size(), order);

  } else if (direction == 1) {

    less_2nd <size_type, value_type> order;
    std::partial_sort(x.begin(), x.begin() + k, x.begin() + x.size(), order);
  }
}

//--------------------------------------------------------------------------------
// Static buffer for partial_argsort, so that we don't have to allocate
// memory each time (faster).
//--------------------------------------------------------------------------------
static SparseVector<size_t, float> partial_argsort_buffer;

//--------------------------------------------------------------------------------
// A partial argsort that can use an already allocated buffer to avoid creating
// a data structure each time it's called. Assumes that the elements to be sorted
// are float, or at least that they have the same size.
//
// A partial sort is much faster than a full sort. The elements after the k first
// in the result are not sorted, except that they are greater (or lesser) than
// all the k first elements. If direction is -1, the sort is in decreasing order.
// If direction is 1, the sort is in increasing order.
//
// The result is returned in the first k positions of the buffer for speed.
//
// Uses a pre-allocated buffer to avoid allocating memory each time a sort
// is needed.
//--------------------------------------------------------------------------------
template<typename InIter, typename OutIter>
inline void partial_argsort(size_t k, InIter begin, InIter end,
                            OutIter sorted, OutIter sorted_end,
                            int direction = -1) {
  {
    assert(0 < k);
    assert(0 < end - begin);
    assert(k <= (size_t) (end - begin));
    assert(k <= (size_t) (sorted_end - sorted));
    assert(direction == -1 || direction == 1);
  }

  typedef size_t size_type;
  typedef float value_type;

  SparseVector <size_type, value_type>& buff = partial_argsort_buffer;

  size_type n = (size_type) (end - begin);

  // Need to clean up, lest the next sort, with a possibly smaller range,
  // picks up values that are not in the current [begin,end).
  buff.resize(n);
  //buff.size() = n;

  InIter it = begin;

  for (size_type i = 0; i != n; ++i, ++it) {
    buff[i].first = i;
    buff[i].second = *it;
  }

  partial_argsort(k, buff, direction);

  for (size_type i = 0; i != k; ++i)
    sorted[i] = buff[i].first;
}

//--------------------------------------------------------------------------------
// QUANTIZE
//--------------------------------------------------------------------------------
template<typename It1>
inline void
update_with_indices_of_non_zeros(unsigned int segment_size,
                                 It1 input_begin, It1 input_end,
                                 It1 prev_begin, It1 prev_end,
                                 It1 curr_begin, It1 curr_end) {
  typedef unsigned int size_type;
  //typedef float value_type;

  size_type input_size = (size_type) (input_end - input_begin);

  std::fill(curr_begin, curr_end, 0);

  for (size_type i = 0; i != input_size; ++i) {

    if (*(input_begin + i) == 0)
      continue;

    size_type begin = i * segment_size;
    size_type end = begin + segment_size;
    bool all_zero = true;

    for (size_type j = begin; j != end; ++j) {

      if (*(prev_begin + j) > 0) {
        all_zero = false;
        *(curr_begin + j) = 1;
      }
    }

    if (all_zero)
      std::fill(curr_begin + begin, curr_begin + end, 1);
  }
}

//--------------------------------------------------------------------------------
// Winner takes all
//--------------------------------------------------------------------------------
/**
 * Finds the maximum in each interval defined by the boundaries, replaces that
 * maximum by a 1, and sets all the other values to 0. Returns the max value
 * over all intervals, and its position.
 */
template<typename I, typename InIter, typename OutIter>
inline void
winnerTakesAll(const std::vector<I>& boundaries, InIter begin1, OutIter begin2) {
  typedef typename std::iterator_traits<InIter>::value_type value_type;

  I max_i = 0, size = (I) boundaries.size();
  value_type max_v = 0;

  for (I i = 0, k = 0; i < size; ++i) {
    max_v = 0;
    max_i = i == 0 ? 0 : boundaries[i - 1];
    while (k < boundaries[i]) {
      if (*begin1 > max_v) {
        max_i = k;
        max_v = *begin1;
      }
      ++k;
      ++begin1;
    }
    *begin2 = (value_type) max_i;
    ++begin2;
  }
}

//--------------------------------------------------------------------------------
/**
 * Winner takes all 2.
 */
template<typename I, typename InIter, typename OutIter>
std::pair<I, typename std::iterator_traits<InIter>::value_type>
winnerTakesAll2(const std::vector<I>& boundaries, InIter begin1, OutIter begin2) {
  I max_i = 0;
  typedef typename std::iterator_traits<InIter>::value_type value_type;
  value_type max_v = 0;

  for (I i = 0, k = 0; i < boundaries.size(); ++i) {
    max_v = 0;
    max_i = i == 0 ? 0 : boundaries[i - 1];
    while (k < boundaries[i]) {
      if (begin1[k] > max_v) {
        begin2[max_i] = 0;
        max_i = k;
        max_v = (value_type) (begin1[k]);
      } else {
        begin2[k] = 0;
      }
      ++k;
    }
    begin2[max_i] = 1;
  }
  return std::make_pair(max_i, max_v);
}

//--------------------------------------------------------------------------------
/**
 * Keeps the values of k winners per segment, where each segment in [begin..end)
 * has length seg_size, and zeroes-out all the other elements.
 * Returns the indices and the values of the winners.
 * For zero segments, we randomly pick a winner, and output its index, with the
 * value zero.
 * If a segment has only zeros, randomly picks a winner.
 */
template<typename I, typename InIter, typename OutIter1, typename OutIter2,
typename RNG>
inline void
winnerTakesAll3(I k, I seg_size, InIter begin, InIter end,
                OutIter1 ind, OutIter2 nz, RNG& rng) {
  typedef I size_type;
  typedef typename std::iterator_traits<InIter>::value_type value_type;

  { // Pre-conditions
    assert(k > 0);
    assert(seg_size > 0);
    assert(k <= seg_size);
    assert((size_type) (end - begin) % seg_size == 0);
  } // End pre-conditions

  typedef select2nd <std::pair<size_type, value_type>> sel2nd;

  InIter seg_begin = begin;
  size_type offset = (size_type) 0;

  for (; seg_begin != end; seg_begin += seg_size, offset += seg_size) {

    InIter seg_end = seg_begin + seg_size;
    size_type offset = (size_type) (seg_begin - begin);

    if (nearlyZeroRange(seg_begin, seg_end)) {

      std::vector<size_type> indices(seg_size);
      random_perm_interval(indices, offset, 1, rng);

      sort(indices.begin(), indices.begin() + k, std::less<size_type>());

      for (size_type i = 0; i != k; ++i, ++ind, ++nz) {
        *ind = indices[i];
        *nz = (value_type) 0;
      }

    } else {

      partial_sort(k, seg_begin, seg_end, ind, nz,
                   predicate_compose<std::greater<value_type>, sel2nd>(),
                   offset, true);
    }
  }
}

////--------------------------------------------------------------------------------
//template <typename I, typename InIter, typename OutIter1, typename OutIter2>
//inline void
//winnerTakesAll3(I k, I seg_size, InIter begin, InIter end,
//                OutIter1 ind, OutIter2 nz)
//{
//  Random rng;
//  winnerTakesAll3(k, seg_size, begin, end, ind, nz, rng);
//}

//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
// IO CONTROL AND MANIPULATORS
//--------------------------------------------------------------------------------
typedef enum {
  CSR = 0, CSR_01, BINARY, AS_DENSE
} SPARSE_IO_TYPE;

struct IOControl {
  int abbr;                  // shorten long vectors output
  bool output_n_elts;        // output vector size at beginning

  bool pair_paren;           // put parens around pairs in vector of pairs
  const char *pair_sep;      // put separator between pair.first and pair.second

  int convert_to_sparse;     // convert dense vector to pos. of non-zeros
  int convert_from_sparse;   // convert from pos. of non-zero to dense 0/1 vector

  SPARSE_IO_TYPE sparse_io;  // do sparse io according to SPARSE_IO_TYPE

  bool bit_vector;           // output 0/1 vector compactly

  inline IOControl(int a = -1, bool s = true, bool pp = false, const char *psep = " ",
                   SPARSE_IO_TYPE smio = CSR,
                   bool cts = false,
                   bool cfs = false,
                   bool bv = false)
  : abbr(a),
    output_n_elts(s),
    pair_paren(pp),
    pair_sep(psep),
    convert_to_sparse(cts),
    convert_from_sparse(cfs),
    sparse_io(smio),
    bit_vector(bv) { }

  inline void reset() {
    abbr = -1;
    output_n_elts = true;
    pair_paren = false;
    pair_sep = " ";
    convert_to_sparse = false;
    convert_from_sparse = false;
    sparse_io = CSR;
    bit_vector = false;
  }
};

extern IOControl io_control;

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator,(std::basic_ostream<CharT, Traits>& out_stream, const T& a) {
  return out_stream << ' ' << a;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator,(std::basic_istream<CharT, Traits>& in_stream, T& a) {
  return in_stream >> a;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
operator,(std::basic_ostream<CharT, Traits>& out_stream,
          std::basic_ostream<CharT, Traits>& (*pf)(std::basic_ostream<CharT, Traits>&)) {
  pf(out_stream);
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
p_paren(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.pair_paren = true;
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
psep_comma(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.pair_sep = ",";
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
psep_dot(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.pair_sep = ".";
  return out_stream;
}

struct abbr {
  int n;
  inline abbr(int _n) : n(_n) { }
};

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits>& out_stream, abbr s) {
  io_control.abbr = s.n;
  return out_stream;
}

struct debug {
  int n;
  inline debug(int _n = -1) : n(_n) { }
};

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits>& out_stream, debug d) {
  io_control.abbr = d.n;
  io_control.output_n_elts = false;
  io_control.pair_sep = ",";
  io_control.pair_paren = true;
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_istream<CharT, Traits>&
from_csr_01(std::basic_istream<CharT, Traits>& in_stream) {
  io_control.convert_from_sparse = true;
  return in_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
to_csr_01(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.convert_to_sparse = true;
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_istream<CharT, Traits>&
bit_vector(std::basic_istream<CharT, Traits>& in_stream) {
  io_control.bit_vector = true;
  return in_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
bit_vector(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.bit_vector = true;
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_istream<CharT, Traits>&
general_vector(std::basic_istream<CharT, Traits>& in_stream) {
  io_control.bit_vector = false;
  io_control.convert_to_sparse = false;
  return in_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
general_vector(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.bit_vector = false;
  io_control.convert_to_sparse = false;
  return out_stream;
}

//--------------------------------------------------------------------------------
// SM IO CONTROL
//--------------------------------------------------------------------------------
struct sparse_format_class {
  SPARSE_IO_TYPE format;

  inline sparse_format_class(SPARSE_IO_TYPE f) : format(f) { }
};

inline sparse_format_class
sparse_format(SPARSE_IO_TYPE f) { return sparse_format_class(f); }

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits>& out_stream, sparse_format_class s) {
  io_control.sparse_io = s.format;
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits>& in_stream, sparse_format_class s) {
  io_control.sparse_io = s.format;
  return in_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
as_dense(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.sparse_io = AS_DENSE;
  return out_stream;
}

template<typename CharT, typename Traits>
inline std::basic_ostream<CharT, Traits>&
as_binary(std::basic_ostream<CharT, Traits>& out_stream) {
  io_control.sparse_io = BINARY;
  return out_stream;
}

//--------------------------------------------------------------------------------
// CHECKERS
//--------------------------------------------------------------------------------
template<typename T1>
struct is_positive_checker {
  T1& var;

  inline is_positive_checker(T1& v) : var(v) { }

  template<typename CharT, typename Traits>
  inline void do_check(std::basic_istream<CharT, Traits>& in_stream) {
    double value = 0;
    in_stream >> value;
    if (value < 0) {
      std::cout << "Value out of range: " << value
      << " - Expected positive or zero value"
      << std::endl;
      exit(-1);
    }
    var = (T1) value;
  }
};

template<typename CharT, typename Traits, typename T1>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits>& in_stream, is_positive_checker<T1> cp) {
  cp.do_check(in_stream);
  return in_stream;
}

template<typename T1>
inline is_positive_checker<T1> assert_positive(T1& var) {
  return is_positive_checker<T1>(var);
}

//--------------------------------------------------------------------------------
// BINARY PERSISTENCE
//--------------------------------------------------------------------------------
template<typename It>
inline void binary_save(std::ostream& out_stream, It begin, It end) {
  typedef typename std::iterator_traits<It>::value_type value_type;
  size_t size = (size_t) (end - begin);
  if (size > 0) {
    char *ptr = (char *) &*begin;
    out_stream.write(ptr, (std::streamsize) size * sizeof(value_type));
  }
}

//--------------------------------------------------------------------------------
template<typename It>
inline void binary_load(std::istream& in_stream, It begin, It end) {
  typedef typename std::iterator_traits<It>::value_type value_type;
  size_t size = (size_t) (end - begin);
  if (size > 0) {
    char *ptr = (char *) &*begin;
    in_stream.read(ptr, (std::streamsize) size * sizeof(value_type));
  }
}

//--------------------------------------------------------------------------------
template<typename T>
inline void binary_save(std::ostream& out_stream, const std::vector<T>& v) {
  binary_save(out_stream, v.begin(), v.end());
}

//--------------------------------------------------------------------------------
template<typename T>
inline void binary_load(std::istream& in_stream, std::vector<T>& v) {
  binary_load(in_stream, v.begin(), v.end());
}

//--------------------------------------------------------------------------------
// STL STREAMING OPERATORS
//--------------------------------------------------------------------------------

//--------------------------------------------------------------------------------
// std::pair
//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& out_stream, const std::pair<T1, T2>& p) {
  if (io_control.pair_paren)
    out_stream << "(";
  out_stream << p.first;
  out_stream << io_control.pair_sep;
  out_stream << p.second;
  if (io_control.pair_paren)
    out_stream << ")";
  return out_stream;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline std::istream& operator>>(std::istream& in_stream, std::pair<T1, T2>& p) {
  in_stream >> p.first >> p.second;
  return in_stream;
}

//--------------------------------------------------------------------------------
// std::vector
//--------------------------------------------------------------------------------
template<typename T, bool>
struct vector_loader {
  inline void load(size_t, std::istream&, std::vector<T>&);
};

//--------------------------------------------------------------------------------
/**
 * Partial specialization of above functor for primitive types.
 */
template<typename T>
struct vector_loader<T, true> {
  inline void load(size_t n, std::istream& in_stream, std::vector<T>& v) {
    if (io_control.convert_from_sparse == CSR_01) {

      std::fill(v.begin(), v.end(), (T) 0);

      for (size_t i = 0; i != n; ++i) {
        int index = 0;
        in_stream >> index;
        v[index] = (T) 1;
      }

    } else if (io_control.bit_vector) {

      for (size_t i = 0; i != n; ++i) {
        float x = 0;
        in_stream >> x;
        if (x)
          v[i] = 1;
        else
          v[i] = 0;
      }

    } else {
      for (size_t i = 0; i != n; ++i)
        in_stream >> v[i];
    }
  }
};

//--------------------------------------------------------------------------------
/**
 * Partial specialization for non-primitive types.
 */
template<typename T>
struct vector_loader<T, false> {
  inline void load(size_t n, std::istream& in_stream, std::vector<T>& v) {
    for (size_t i = 0; i != n; ++i)
      in_stream >> v[i];
  }
};

//--------------------------------------------------------------------------------
/**
 * Factory that will instantiate the right functor to call depending on whether
 * T is a primitive type or not.
 */
//template<typename T>
//inline void vector_load(size_t n, std::istream& in_stream, std::vector<T>& v) {
//  vector_loader <T, boost::is_fundamental<T>::value> loader;
//  loader.load(n, in_stream, v);
//}

//--------------------------------------------------------------------------------
template<typename T, bool>
struct vector_saver {
  inline void save(size_t n, std::ostream& out_stream, const std::vector<T>& v);
};

//--------------------------------------------------------------------------------
/**
 * Partial specialization for primitive types.
 */
template<typename T>
struct vector_saver<T, true> {
  inline void save(size_t n, std::ostream& out_stream, const std::vector<T>& v) {
    if (io_control.output_n_elts)
      out_stream << n << ' ';

    if (io_control.abbr > 0)
      n = std::min((size_t) io_control.abbr, n);

    if (io_control.convert_to_sparse == CSR_01) {

      for (size_t i = 0; i != n; ++i)
        if (!is_zero(v[i]))
          out_stream << i << ' ';

    } else if (io_control.bit_vector) {

      size_t k = 7;
      for (size_t i = 0; i != v.size(); ++i) {
        out_stream << (is_zero(v[i]) ? '0' : '1');
        if (i == k) {
          out_stream << ' ';
          k += 8;
        }
      }

    } else {

      for (size_t i = 0; i != n; ++i)
        out_stream << v[i] << ' ';
    }

    if (io_control.abbr > 0 && n < v.size()) {
      size_t rest = v.size() - n;
      out_stream << "[+" << rest << "/" << count_non_zeros(v) << "]";
    }
  }
};

//--------------------------------------------------------------------------------
/**
 * Partial specialization for non-primitive types.
 */
template<typename T>
struct vector_saver<T, false> {
  inline void save(size_t n, std::ostream& out_stream, const std::vector<T>& v) {
    if (io_control.output_n_elts)
      out_stream << n << ' ';

    if (io_control.abbr > 0)
      n = std::min((size_t) io_control.abbr, n);

    for (size_t i = 0; i != n; ++i)
      out_stream << v[i] << ' ';

    if (io_control.abbr > 0 && n < v.size()) {
      size_t rest = v.size() - n;
      out_stream << "[+" << rest << "/" << count_non_zeros(v) << "]";
    }
  }
};

//--------------------------------------------------------------------------------
/**
 * Factory that will instantiate the right functor to call depending on whether
 * T is a primitive type or not.
 */
//template<typename T>
//inline void vector_save(size_t n, std::ostream& out_stream, const std::vector<T>& v) {
//  vector_saver <T, boost::is_fundamental<T>::value> saver;
//  saver.save(n, out_stream, v);
//}

//--------------------------------------------------------------------------------
/**
 * Saves the size of the vector.
 */
//template<typename T>
//inline std::ostream& operator<<(std::ostream& out_stream, const std::vector<T>& v) {
//  vector_save(v.size(), out_stream, v);
//  return out_stream;
//}

//--------------------------------------------------------------------------------
template<typename T>
inline std::ostream& operator<<(std::ostream& out_stream, const std::list<T>& l) {
  for (typename std::list<T>::const_iterator i = l.begin(); i != l.end(); ++i)
    out_stream << *i << " ";
  return out_stream;
}

//--------------------------------------------------------------------------------
/**
 * Reads in size of the vector, and redimensions it, except if we are reading
 * a sparse binary vector.
 */
template<typename T>
inline std::istream&
operator>>(std::istream& in_stream, std::vector<T>& v) {
  size_t n = 0;
  in_stream >> n;
  v.resize(n);
  vector_load(n, in_stream, v);
  return in_stream;
}

//--------------------------------------------------------------------------------
// std::set
//--------------------------------------------------------------------------------
template<typename T1>
inline std::ostream& operator<<(std::ostream& out_stream, const std::set<T1>& m) {
  typename std::set<T1>::const_iterator
  it = m.begin(), end = m.end();

  while (it != end) {
    out_stream << *it << ' ';
    ++it;
  }

  return out_stream;
}

//--------------------------------------------------------------------------------
// std::map
//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& out_stream, const std::map<T1, T2>& m) {
  out_stream << m.size() << " ";

  typename std::map<T1, T2>::const_iterator
  it = m.begin(), end = m.end();

  while (it != end) {
    out_stream << it->first << ' ' << it->second << ' ';
    ++it;
  }

  return out_stream;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline std::istream& operator>>(std::istream& in_stream, std::map<T1, T2>& m) {
  int size = 0;
  in_stream >> size;

  for (int i = 0; i != size; ++i) {
    T1 k;
    T2 v;
    in_stream >> k >> v;
    m.insert(std::make_pair(k, v));
  }

  return in_stream;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline std::ostream& operator<<(std::ostream& out_stream, const std::multimap<T1, T2>& m) {
  out_stream << m.size() << " ";

  typename std::multimap<T1, T2>::const_iterator
  it = m.begin(), end = m.end();

  while (it != end) {
    out_stream << it->first << ' ' << it->second << ' ';
    ++it;
  }

  return out_stream;
}

//--------------------------------------------------------------------------------
template<typename T1, typename T2>
inline std::istream& operator>>(std::istream& in_stream, std::multimap<T1, T2>& m) {
  int size = 0;
  in_stream >> size;

  for (int i = 0; i != size; ++i) {
    T1 k;
    T2 v;
    in_stream >> k >> v;
    m.insert(std::make_pair(k, v));
  }

  return in_stream;
}

//--------------------------------------------------------------------------------
// MISCELLANEOUS
//--------------------------------------------------------------------------------
template<typename T>
inline void show_all_differences(const std::vector<T>& x, const std::vector<T>& y) {
  std::vector<size_t> diffs;
  find_all_differences(x, y, diffs);
  std::cout << diffs.size() << " differences: " << std::endl;
  for (size_t i = 0; i != diffs.size(); ++i)
    std::cout << "(at:" << diffs[i]
    << " y=" << x[diffs[i]]
    << ", ans=" << y[diffs[i]] << ")";
  std::cout << std::endl;
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
IOControl io_control;

//---------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
template<typename T>
inline void print_int_w_comma(std::ostream& out_stream, const T& x) {
  std::stringstream s;
  s << x;
  const std::string str = s.str();
  std::vector<char> v;
  for (int i = str.size() - 1; 0 <= i; --i) {
    v.push_back(str[i]);
    if (0 < i && i % 3 == 0)
      v.push_back(',');
  }
  std::reverse(v.begin(), v.end());
  for (size_t i = 0; i != v.size(); ++i)
    out_stream << v[i];
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// Issues with boost::binomial_coefficient, and this one will work also for
// real n, k. Also, lgamma and al. are in math.h, no need to have boost.
// But watch out!!! Call with float or double only!!!
// Choose(4, 2) = 5 !! WRONG
// Choose(4.0, 2.0) = 6 CORRECT
// Doesn't work too well for really large numbers, runs into precision problems
// quickly.
template<typename T1, typename T2>
long double Choose(T1 _n, T2 _k) {
  long double n = (long double) _n;
  long double k = (long double) _k;
  return exp2l(lgammal(n + 1.0) - lgammal(k + 1.0) - lgammal(n - k + 1.0));
}

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// requires 64 bits flags, because lib was built for 64 bits...
/*
mpz_class gmp_factorial(size_t n)
{
  if (n == 0)
    return 1;
  mpz_class f = 1;
  for (size_t i = 2; i <= n; ++i)
    f *= i;
  return f;
}

//--------------------------------------------------------------------------------
mpz_class gmp_binomial(size_t n, size_t k)
{
  return gmp_factorial(n) / (gmp_factorial(k) * gmp_factorial(n - k));
}

//--------------------------------------------------------------------------------
mpz_class gmp_multinomial(const std::vector<size_t>& n)
{
  mpz_class den = 1;
  for (size_t i = 0; i != n.size(); ++i)
    den *= gmp_factorial(n[i]);
  return gmp_factorial(sum(n)) / den;
}
*/

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------
// Enumerate multi-index whose each element is < the corresponding bound:
//
// vector<size_t> b;
// b += 3,2,3;
// enumerate e(b);
// while (e.ok())
//   cout << e.next() << endl;
//
// gives:
//
// 0 0 0
// 0 0 1
// 0 0 2
// 0 1 0
// 0 1 1
// 0 1 2
// 1 0 0
// 1 0 1
// 1 0 2
// 1 1 0
// 1 1 1
// 1 1 2
// 2 0 0
// 2 0 1
// 2 0 2
// 2 1 0
// 2 1 1
// 2 1 2
struct enumerate {
  size_t i, n;
  std::vector<size_t> bounds;
  std::vector<size_t> current;

  inline enumerate(const std::vector<size_t>& _bounds)
  : i(0), n(product(_bounds)),
    bounds(_bounds),
    current(bounds.size(), 0) { }

  inline bool ok() { return i < n; }

  inline std::vector<size_t>& next() {
    if (i == 0) {
      ++i;
      return current;
    }
    int j = bounds.size() - 1;
    ++current[j];
    while (current[j] == bounds[j]) {
      current[j] = 0;
      ++current[--j];
    }
    ++i;
    return current;
  }
};

//--------------------------------------------------------------------------------
template<typename T>
inline T gcd(T a, T b) {
  T g;
  while (b != 0) {
    g = b;
    b = a % b;
    a = g;
  }
  return a;
}

//--------------------------------------------------------------------------------
// Integer power
template<typename T1, typename T2>
inline T1 ipow(T1 x, T2 p) {
  T1 y = 1;
  while (p > 0) {
    if (p & 1) {
      y *= x;
      p--;
    } else {
      x *= x;
      p >>= 1;
    }
  }
  return y;
}

//--------------------------------------------------------------------------------
/*
 * A utility function that visualizes a vector of non-negative values.
 * For debugging.
 */
template<typename T>
void PlotNonNegativeVector(const std::vector<T>& x, size_t W = 80) {
  T M = max_value(x);
  if (M <= 0)
    return;

  float ratio = float(W) / float(M);
  size_t w = int(log(x.size()) / log(10)) + 1;

  for (size_t i = 0; i < x.size(); ++i) {
    std::cout << std::setw(w) << i << ":";
    size_t upTo = ratio * x[i];
    for (size_t j = 0; j < upTo; ++j)
      std::cout << "+";
    std::cout << std::endl;
  }
}

//--------------------------------------------------------------------------------
int to_int(unsigned char s[]) {
  int v = 0;
  for (int k = 0; k < 8; ++k)
    v += (s[k] == '1') * (1 << (7 - k));
  return v;
}

//--------------------------------------------------------------------------------
void to_dec(int n, std::string& d) {
  for (int i = 7; i >= 0; --i)
    d[7 - i] = n & (1 << i) ? '1' : '0';
}

//--------------------------------------------------------------------------------
int hamming(int a, int b) {
  int h = 0;
  for (int i = 7; i >= 0; --i)
    h += (a & (1 << i)) == (b & (1 << i));
  return 8 - h;
}

//--------------------------------------------------------------------------------
int hamming(unsigned char a, unsigned char b) {
  int h = 0;
  for (int i = 7; i >= 0; --i)
    h += (a & (1 << i)) == (b & (1 << i));
  return 8 - h;
}


//----------------------------------------------------------------------------------------------------------------------
///**
// * Compile-time power of 2 function.
// */
//inline constexpr uint64_t pow2(uint8_t exponent)
//{
//  return (exponent == 0) ? 1 : (2 * pow2(exponent-uint8_t(1)));
//}

//----------------------------------------------------------------------------------------------------------------------
// Extracts size of char[]
template<typename T, std::size_t N>
inline constexpr std::size_t size(T (&)[N]) noexcept {
  return N;
}

////----------------------------------------------------------------------------------------------------------------------
//inline void dump_memory(void* ptr, std::size_t size, std::ostream& os = std::cout)
//{
//  typedef unsigned char byte;
//  typedef unsigned long uint;
//
//  // Allow direct arithmetic on the pointer
//  uint iptr = reinterpret_cast<uint>(ptr);
//
//  os << "-----------------------------------------------------------------------\n";
//  os << boost::format("%d bytes") % size;
//
//  // Get number of digits
//  uint indent = (uint) std::log10(size) + 1;
//
//  // Write the address offsets along the top row
//  // Account for the indent of "X bytes"
//  os << std::string(13 - indent, ' ');
//  for(std::size_t i = 0; i < 16; ++i)
//  {
//    if(i %  4 == 0){os << " ";}        // Spaces between every 4 bytes
//    os << boost::format(" %2hhX") % i; // Write the address offset
//  }
//
//  // If the object is not aligned
//  if(iptr % 16 != 0)
//  {
//    // Print the first address
//    os << boost::format("\n0x%016lX:") % (iptr & ~15);
//
//    // Indent to the offset
//    for(std::size_t i = 0; i < iptr % 16; ++i)
//    {
//      os << "   ";
//      if(i % 4 == 0){os << " ";}
//    }
//  }
//
//  // Dump the memory
//  for(std::size_t i = 0; i < size; ++i, ++iptr)
//  {
////    if(iptr % 16 == 0) {
////      if (i >= 16) {
////        os << " ";
////        char s[16];
////        std::copy((byte *) iptr - 16, (byte *) iptr, s);
////        os << " " << s;
////      }
////    }
//
//    // New line and address every 16 bytes, spaces every 4 bytes
//    if(iptr % 16 == 0) {
//      os << boost::format("\n0x%016lX:") % iptr;
//    }
//    if(iptr %  4 == 0){os << " ";}
//
//    // Write the address contents
//    os << boost::format(" %02hhX")
//          % static_cast<uint>(*reinterpret_cast<byte*>(iptr));
//  }
//
////  os << " ";
////  char s[16];
////  std::copy((byte *) iptr - 16, (byte *) iptr, s);
////  os << " " << s;
//
//  os << "\n-----------------------------------------------------------------------"
//  << std::endl;
//}
//
////----------------------------------------------------------------------------------------------------------------------
//inline void dump_memory_with_context(void* ptr,
//                                     std::size_t size,
//                                     std::ostream& os = std::cout)
//{
//  // Allow direct arithmetic on the pointer
//  unsigned long sptr = reinterpret_cast<unsigned long>(ptr); // Start pointer
//  unsigned long eptr = sptr + size;                          // End pointer
//
//  sptr &= ~15; // Round down to the last multiple of 16
//  sptr -=  16; // Step back one line for context
//  eptr &= ~15; // Round down to the last multiple of 16
//  eptr +=  32; // Step forward one line for context
//
//  // Dump memory
//  dump_memory(reinterpret_cast<void*>(sptr), eptr - sptr, os);
//}

////----------------------------------------------------------------------------------------------------------------------
//inline void TEST(bool expr, const char* msg =nullptr) {
//  if (!expr) {
//    void* callstack[128];
//    int frames = backtrace(callstack, 128);
//    cerr << "Test failed";
//    if (msg)
//      cerr << msg;
//    cerr << endl;
//    char** bt = backtrace_symbols(callstack, frames);
//    for (int i = 0; i < frames; ++i)
//      cerr << bt[frames-i-1] << endl;
//    exit(-1);
//  }
//}

//--------------------------------------------------------------------------------
} // end namespace utils
#endif //CPP_UTILS_TOOLBOX_HPP
