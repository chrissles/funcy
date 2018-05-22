/* MLPACK 0.2
 *
 * Copyright (c) 2008, 2009 Alexander Gray,
 *                          Garry Boyer,
 *                          Ryan Riegel,
 *                          Nikolaos Vasiloglou,
 *                          Dongryeol Lee,
 *                          Chip Mappus, 
 *                          Nishant Mehta,
 *                          Hua Ouyang,
 *                          Parikshit Ram,
 *                          Long Tran,
 *                          Wee Chin Wong
 *
 * Copyright (c) 2008, 2009 Georgia Institute of Technology
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */
/**
 * @file intmap.h
 *
 * Dense integer-to-value map.
 */

#ifndef COL_INTMAP_H
#define COL_INTMAP_H

#include "../base/base.h"

/**
 * A dense grow-as-needed array that serves as an integer-keyed map.
 */
template<class TValue>
class DenseIntMap {
 public:
  typedef TValue Value;

 private:
  Value *ptr_;
  fl__index_t size_;
  Value default_value_;

  OT_DEF(DenseIntMap) {
    OT_MY_OBJECT(size_);
    OT_MY_OBJECT(default_value_);
    OT_MALLOC_ARRAY(ptr_, size_);
  }

 public:
  /** Creates a blank mapping. */
  void Init() {
    ptr_ = NULL;
    size_ = 0;
  }

  /** Accesses the default value.  Use this to set it. */
  Value& default_value() {
    return default_value_;
  }
  /** Accesses the default value. */
  const Value& default_value() const {
    return default_value_;
  }

  /**
   * Gets a non-inclusive upper bound on the last non-default element.
   *
   * Iterating up to size() is guaranteed to hit all elements without
   * growing the internal array.
   */
  fl__index_t size() const {
    return size_;
  }

  /**
   * Accesses an element, expanding if necessary.
   *
   * If you are just probing, beware that this might actually grow the array!
   */
  Value& operator [] (fl__index_t index) {
    DEBUG_BOUNDS(index, BIG_BAD_NUMBER);
    if ((index >= size_)) {
      fl__index_t old_size = size_;
      size_ = std::max(size_ * 2, index + 1);
      ptr_ = mem::Realloc(ptr_, size_);
      for (fl__index_t i = old_size; i < size_; i++) {
        new(ptr_+i)Value(default_value_);
      }
    }
    return ptr_[index];
  }
  /**
   * Accesses an element from a static context.
   */
  const Value& operator [] (fl__index_t index) const {
    return get(index);
  }
  /**
   * Accesses an element, never growing the internal representation.
   */
  const Value& get(fl__index_t index) const {
    DEBUG_BOUNDS(index, BIG_BAD_NUMBER);
    if ((index < size_)) {
      return ptr_[index];
    } else {
      return default_value_;
    }
  }
};

#endif
