//
// This file is part of the CAMORRA library.
// Copyright (C) 2010 Gijs van den Oord.
// CAMORRA is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PPRANKNUM_H_
#define CAMGEN_PPRANKNUM_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Some useful definitions of function calls that will make coding easier... *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define RANK_NUMERAL_BASE0	first
#define RANK_NUMERAL_BASE1	second
#define RANK_NUMERAL_BASE2	third
#define RANK_NUMERAL_BASE3	fourth

#define RANK_NUMERAL_LVL1(i) RANK_NUMERAL_BASE##i

#define RANK_NUMERAL(i) RANK_NUMERAL_LVL1(i)

#endif /*CAMGEN_PPRANKNUM_H_*/

