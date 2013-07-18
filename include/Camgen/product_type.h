//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PRODUCT_TYPE_H_
#define CAMGEN_PRODUCT_TYPE_H_

#include <complex>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Class template selecting the numerical type of a product of real and/or *
 * complex number:                                                         *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
	template<class value_t,class type1,class type2>class product_type{};
	
	template<class value_t>class product_type<value_t,value_t,value_t>
	{
		public:
			typedef value_t value_type;
	};
	
	template<class value_t>class product_type<value_t,std::complex<value_t>,value_t>
	{
		public:
			typedef std::complex<value_t> value_type;
	};

	template<class value_t>class product_type<value_t,value_t,std::complex<value_t> >
	{
		public:
			typedef std::complex<value_t> value_type;
	};
	
	template<class value_t>class product_type<value_t,std::complex<value_t>,std::complex<value_t> >
	{
		public:
			typedef std::complex<value_t> value_type;
	};
}

#endif /*CAMGEN_PRODUCT_TYPE_H_*/

