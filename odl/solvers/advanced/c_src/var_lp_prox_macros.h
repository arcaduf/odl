/* Copyright 2014-2016 The ODL development group

 This file is part of ODL.

 ODL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ODL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ODL. If not, see <http://www.gnu.org/licenses/>.
*/

// The SINGLE or DOUBLE macros must be defined externally

#include <math.h>

#if defined SINGLE
#define FLOAT float
#define MAX fmaxf
#define MIN fminf
#define ABS fabsf
#define POW powf
#define IMPL(fname) fname ## __float

#elif defined DOUBLE
#define FLOAT double
#define MAX fmax
#define MIN fmin
#define ABS fabs
#define POW pow
#define IMPL(fname) fname ## __double

#else
#error "Either SINGLE or DOUBLE must be defined!"

#endif

