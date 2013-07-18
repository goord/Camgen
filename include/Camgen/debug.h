//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DEBUG_H_
#define CAMGEN_DEBUG_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Preprocessor macro definitions for debug logging in Camgen.   *
 *                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <iostream>
#include <cstdlib>

/* Uniform debug flag setup: */

#if defined DEBUG
#define CAMGEN_DEBUG
#elif defined _DEBUG
#define CAMGEN_DEBUG
#endif

#ifdef NDEBUG
#undef CAMGEN_DEBUG
#endif

/* Selecting compiler's function name macro, if existing: */

#if defined(__PRETTY_FUNCTION__)
#define __CAMGEN_FUNC__ __PRETTY_FUNCTION__
#elif defined(__FUNCTION__)
#define __CAMGEN_FUNC__ __FUNCTION__
#elif defined(__func__)
#define __CAMGEN_FUNC__ __func__
#elif defined(_MSC_VER)
#define __CAMGEN_FUNC__ __FUNCTION__
#elif defined(__GNUC__)
#define __CAMGEN_FUNC__ __PRETTY_FUNCTION__
#elif defined(__BORLANDC__)
#define __CAMGEN_FUNC__ __FUNC__
#else
#warning "no function name printing facility detected...proceeding with <unknown> function names"
#define __CAMGEN_FUNC__ "<unknown>"
#endif

/* Helper macro definitions: */

#define CAMGEN_STREAMLOC "file "<<__FILE__<<", function "<<__CAMGEN_FUNC__<<", line "<<__LINE__<<": "

/* Message macro definition: */

#ifdef CAMGEN_DEBUG
#define CAMGEN_MESSAGE(msg) std::cerr<<CAMGEN_STREAMLOC<<msg<<std::endl
#else
#define CAMGEN_MESSAGE(msg) (void)0
#endif

/* Conditional message macro definition: */

#ifdef CAMGEN_DEBUG
#define CAMGEN_MESSAGE_IF(cond,msg) 		\
do						\
{						\
    if(cond)					\
    {						\
	CAMGEN_MESSAGE(msg);			\
    }						\
}						\
while(0)
#else
#define CAMGEN_MESSAGE_IF(cond,msg) (void)0
#endif

/* Error macro definition: */

#ifdef CAMGEN_DEBUG
#define CAMGEN_ERROR(msg) 			\
do						\
{						\
    CAMGEN_MESSAGE(msg);			\
    std::abort();				\
}						\
while(0)
#else
#define CAMGEN_ERROR(msg) (void)0
#endif

/* Conditional error macro definition: */

#ifdef CAMGEN_DEBUG
#define CAMGEN_ERROR_IF(cond,msg) 		\
do						\
{						\
    if(cond)					\
    {						\
	CAMGEN_MESSAGE(msg);			\
	std::abort();				\
    }						\
}						\
while(0)
#else
#define CAMGEN_ERROR_IF(cond,msg) (void)0
#endif

#endif /*CAMGEN_DEBUG_H_*/

