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
 * @file base.h
 *
 * Includes all of FASTlib's bare necessities.
 *
 * @see common.h
 * @see debug.h
 * @see compiler.h
 * @see cc.h
 * @see otrav.h
 */

#ifndef BASE_H
#define BASE_H



//#include "common.h"
//#include "debug.h"
#include "common.h"
#include "debug.h"

#ifdef __cplusplus
#include "cc.h"
#include "ccmem.h"
#include "otrav.h"
//#include "cc.h"
//#include "ccmem.h"
//#include "otrav.h"
#endif

#if 0
#ifdef __GNUC__
#ifdef __MINGW32__
#include <float.h>
#define isnan(x) _isnan(x)
#endif
#endif
#endif

#endif /* BASE_H */
