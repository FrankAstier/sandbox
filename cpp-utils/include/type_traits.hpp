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

#ifndef CPP_UTILS_TYPE_TRAITS_HPP
#define CPP_UTILS_TYPE_TRAITS_HPP

#include <type_traits>
#include <vector>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>

/**
 * Additional type traits that can be useful.
 */
namespace utils {

template<typename T>
struct is_std_container : public std::false_type {};

template<typename T, class Alloc>
struct is_std_container<std::vector<T, Alloc>> : public std::true_type {};

template<typename T, class Alloc>
struct is_std_container<std::list<T, Alloc>> : public std::true_type {};

template<typename K, typename V, class Compare, class Alloc>
struct is_std_container<std::map<K, V, Compare, Alloc>> : public std::true_type {};

template<typename K, class Compare, class Alloc>
struct is_std_container<std::set<K, Compare, Alloc>> : public std::true_type {};

template<typename K, typename V, class Hash, class Pred, class Alloc>
struct is_std_container<std::unordered_map<K, V, Hash, Pred, Alloc>> : public std::true_type {};

template<typename K, class Hash, class Pred, class Alloc>
struct is_std_container<std::unordered_set<K, Hash, Pred, Alloc>> : public std::true_type {};

} // end namespace utils

#endif //CPP_UTILS_TYPE_TRAITS_HPP
