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
 * @file deprecated.h
 *
 * Definitions and wrappers so outdated symbols may still be used.
 */

#ifndef BASE_DEPRECATED_H
#define BASE_DEPRECATED_H

#include "base.h"

/* common.h */
#define print_notify_headers print_notify_locs
void percent_indicator(const char *name, uint64 num, uint64 den) {
  fl_print_progress(name, num * 100 / den);
}

#define SUCCESS_FROM_INT(x) SUCCESS_FROM_C(x)

/* compiler.h */
#define EXTERN_C_START EXTERN_C_BEGIN
#define IS_CONSTANT_EXPRESION(expr) IS_CONST_EXPR(expr)
#define COMPILER_NORETURN COMPILER_NO_RETURN
#define COMPILER_NOINLINE COMPILER_NO_INLINE

/* debug.h */
#define debug_verbosity verbosity_level
#define DEBUG_GOT_HERE(min_verbosity) \
    VERBOSE_GOT_HERE(min_verbosity)

#ifdef __cplusplus
/* cc.h */
#define FORBID_COPY(C) FORBID_ACCIDENTAL_COPIES(C)
#define CC_ASSIGNMENT_OPERATOR(C) ASSIGN_VIA_COPY_CONSTRUCTION(C)
#define DEFINE_INEQUALITY_COMPARATORS(C) EXPAND_LESS_THAN(C)
#define DEFINE_ALL_COMPARATORS(C) EXPAND_LESS_THAN(C) EXPAND_EQUALS(C)
#define DEFINE_INEQUALITY_COMPARATORS_HETERO(C, T) \
    EXPAND_HETERO_LESS_THAN(C, T) EXPAND_HETERO_LESS_THAN(T, C)
#define DEFINE_ALL_COMPARATORS_HETERO(C, T) \
    DEFINE_INEQUALITY_COMPARATORS_HETERO(C, T) EXPAND_HETERO_EQUALS(C, T)

/* ccmem.h */
namespace mem {
  template<typename T>
  inline T *ZeroBytes(T *array, size_t bytes) {
    return BitZeroBytes(array, bytes);
  }
  template<typename T>
  inline T *Zero(T *array, size_t elems = 1) {
    return BitZero(array, elems);
  }
  template<typename T>
  inline T *AllocZeroed(size_t elems = 1) {
    return AllocBitZero<T>(elems);
  }

  template<typename T, typename U>
  inline T *CopyBytes(T *dest, const U *src, size_t bytes) {
    return BitCopyBytes(dest, reinterpret_cast<const T *>(src), bytes);
  }
  template<typename T, typename U>
  inline T *Copy(T *dest, const U *src, size_t elems = 1) {
    return BitCopy(dest, reinterpret_cast<const T *>(src), elems);
  }
  template<typename T>
  inline T *DupBytes(const T *src, size_t bytes) {
    return AllocBitCopyBytes(src, bytes);
  }
  template<typename T>
  inline T *Dup(const T *src, size_t elems = 1) {
    return AllocBitCopy(src, elems);
  }

  template<typename T>
  inline T *Resize(T *array, size_t elems = 1) {
    return Realloc(array, elems);
  }

  template<typename T>
  inline T *SwapBytes(T *a, T *b, size_t bytes) {
    return BitSwapBytes(a, b, bytes);
  }
  template<typename T>
  inline T *Swap(T *a, T *b, size_t elems = 1) {
    return BitSwap(a, b, elems);
  }

  template<typename T>
  inline T *ConstructAll(T *array, size_t elems) {
    return Construct(array, elems);
  }
  template<typename T, typename U>
  inline T *Construct(T *ptr, const U &init) {
    return RepeatConstruct(ptr, init, 1);
  }
  template<typename T, typename U>
  inline T *ConstructAll(T *array, const U &init, size_t elems) {
    return RepeatConstruct(array, init, elems);
  }
  template<typename T>
  inline T *DestructAll(T *array, size_t elems) {
    return Destruct(array, elems);
  }

  template<typename T, typename U>
  inline T *AllocConstruct(const U &init, size_t elems) {
    return AllocRepeatConstruct<T>(init, elems);
  }
  template<typename T>
  inline T *DupConstruct(const T *src, size_t elems = 1) {
    return AllocCopyConstruct(src, elems);
  }

  template<typename T>
  inline T *PointerAdd(const T *ptr, ptrdiff_t bytes) {
    return PtrAddBytes(ptr, bytes);
  }
  template<typename T, typename U>
  inline ptrdiff_t *PointerDiff(const T *lhs, const U *rhs) {
    return PtrDiffBytes(lhs, rhs);
  }
  template<typename T>
  inline ptrdiff_t *PointerAbsoluteAddress(const T *ptr) {
    return PtrAbsAddr(ptr);
  }
  template<typename T, typename U>
  inline bool PointersEqual(const T *lhs, const U *rhs) {
    return PtrsEqual(lhs, rhs);
  }
};
#endif

#endif /* BASE_DEPRECATED_H */
