//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_DEF_ARGS_H_
#define CAMGEN_DEF_ARGS_H_

#include <Camgen/type_holders.h>

/* Macro's facilitating the passing of type definitions: */

#define DEFINE_MODEL_TYPE(M)	typedef M model_type

#define DEFINE_SIZE_TYPE(M)												\
typedef typename get_basic_types<M>::size_type size_type

#define DEFINE_VALUE_TYPE(M)												\
typedef typename get_basic_types<M>::value_type value_type

#define DEFINE_R_VALUE_TYPE(M)												\
typedef typename get_basic_types<M>::r_value_type r_value_type

#define DEFINE_TENSOR_TYPE(M)												\
typedef typename get_basic_types<M>::tensor_type tensor_type

#define DEFINE_ITERATOR(M)												\
typedef typename get_basic_types<M>::iterator iterator

#define DEFINE_VECTOR_TYPE(M)												\
typedef typename get_basic_types<M>::vector_type vector_type

#define DEFINE_MOMENTUM_TYPE(M)												\
typedef typename get_basic_types<M>::momentum_type momentum_type

#define DEFINE_SPACETIME_TYPE(M)											\
typedef typename get_spacetime_type<M,M::dimension>::type spacetime_type

#define DEFINE_DIRAC_ALGEBRA_TYPE(M)											\
typedef typename get_Dirac_algebra_type<M>::type Dirac_algebra_type

#define DEFINE_GROUP_TYPE(M,G)												\
typedef typename get_implementation<M,G>::type group_type

#define DEFINE_REP_TYPE(M,R)												\
typedef typename get_implementation<M,R>::type rep_type

#define DEFINE_BASIC_TYPES(M)												\
typedef typename get_basic_types<M>::model_type model_type;								\
typedef typename get_basic_types<M>::size_type size_type;								\
typedef typename get_basic_types<M>::r_value_type r_value_type;								\
typedef typename get_basic_types<M>::value_type value_type;								\
typedef typename get_basic_types<M>::tensor_type tensor_type;								\
typedef typename get_basic_types<M>::iterator iterator;									\
typedef typename get_basic_types<M>::const_iterator const_iterator;							\
typedef typename get_basic_types<M>::momentum_type momentum_type

/* Macro defining the argument list of vertex class methods: */

#define ARG_LIST const value_type& factor,const std::vector<const value_type*>& couplings,std::vector<iterator>& iters,const std::vector<const momentum_type*>& momenta

/* Macro defining the argument list of vertex class methods in the case where
 * the program keeps track of nonzero propagating colour modes: */

#define CFD_ARG_LIST const value_type& factor,const std::vector<const value_type*>& couplings,std::vector<iterator>& iters,const std::vector<const momentum_type*>& momenta,std::set<iterator>& produced_iters

/* Preprocessor definitions of the arguments in the lists above: */

#define C_0 (*couplings[0])
#define C_1 (*couplings[1])
#define C_2 (*couplings[2])
#define A_0 iters[0]
#define A_1 iters[1]
#define A_2 iters[2]
#define A_3 iters[3]
#define P_0 (*momenta[0])
#define P_1 (*momenta[1])
#define P_2 (*momenta[2])
#define P_3 (*momenta[3])

#endif /*CAMGEN_DEF_ARGS_H_*/

