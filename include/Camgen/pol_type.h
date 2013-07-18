//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_POL_TYPE_H_
#define CAMGEN_POL_TYPE_H_

#include <Camgen/vector.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the polarised spin vector class, needed for the *
 * construction of massive spinors. The choice of this spin vector is defined as *
 *                                                                               *
 * 			s = e - (e.p/(m^2+m*E))*(p + p0)                         *
 *                                                                               *
 * where e is the spacelike vector in the polarisation direction and p0 the      *
 * rest-frame momentum. Note that in the rest-frame, p = p0, we have s = e.      *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Wrapper class for the spin vector type: */
    
    class polarised_type
    {
	public:
	    
	    /* Implementation class template, taking the spacetem type and beam
	     * direction as template arguments: */
	    
	    template<class spacetime_t,int beam_dir>class implementation
	    {
		public:
		    typedef vector<typename spacetime_t::value_type,spacetime_t::base_type::index_range> vector_type;

		    /* Construction of the spin vector. Simply calls a static
		     * member fuction of the spacetime implementation class
		     * template. Takes as arguments a vector reference to fill,
		     * a pointer to the momentum and a pointer to the mass: */
		    
		    static vector_type& make_spinvec(vector_type& s,const vector_type* p,const typename spacetime_t::value_type* m)
		    {
			return spacetime_t::template make_polarised_spinvec<beam_dir>(s,*p,*m);
		    }
	    };
    };
    
/* If the massive_spinor_factory class has been declared already,
 * appropriate specialisations are defined, for he purpose of optmisation: */
}

#ifdef CAMGEN_M_SPINOR_FAC_H_

#ifdef CAMGEN_PAULI_BASIS_H_
#include <Camgen/polspinPb.h>
#endif /*CAMGEN_PAULI_BASIS_H_*/

#ifdef CAMGEN_WEYL_BASIS_H_
#include <Camgen/polspinWb.h>
#endif /*CAMGEN_WEYL_BASIS_H_*/

#endif /*CAMGEN_M_SPINOR_FAC_H_*/

#endif /*CAMGEN_POL_TYPE_H_*/


