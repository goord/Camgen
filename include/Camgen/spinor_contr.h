//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_SPINOR_CONTR_H_
#define CAMGEN_SPINOR_CONTR_H_

#include <Camgen/eval.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the spinor contraction class for the evaluation *
 * of the final-current contraction if this current is fermionic. The evaluation *
 * specialisation contains the extra functions Cc_first, etc. to reverse the     *
 * first current if a clash occurs between fermion flows (which can happen if    *
 * the contracted wave function belongs to an incoming Majorana fermion).        *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class spinor_contraction
    {
	private:

	    /* Dirac algebra type definition: */

	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	public:

	    /* Fermionic property flag: */

	    static const bool fermionic=true;
	    
	    /* Contraction tensor size: */
	    
	    static const std::size_t tensor_size=(Dirac_algebra_type::base_type::index_range)*(Dirac_algebra_type::base_type::index_range);
	    
	    /* Contracted subamplitude size: */
	    
	    static const std::size_t size=Dirac_algebra_type::base_type::index_range;
    };
    template<class model_t>const bool spinor_contraction<model_t>::fermionic;
    template<class model_t>const std::size_t spinor_contraction<model_t>::tensor_size;

    /* Specialisation of the evaluate class template for spinor contractions: */

    template<class model_t>class evaluate< spinor_contraction<model_t> >
    {
	public:
	    
	    /* The usual type definitions for evaluate class specialisations: */

	    DEFINE_BASIC_TYPES(model_t);

	private:	
	    
	    /* Necessary Dirac algebra type definition: */
	    
	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	public:
	    
	    /* Initialisation phase: */

	    static void initialise()
	    {
		Dirac_algebra_type::initialise();
	    }
	    
	    /* Index range vector filling function for checking purposes: */
	    
	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		r.resize(1);
		r[0]=Dirac_algebra_type::base_type::index_range;
		return r;
	    }

	    /* Contraction algorithm in the discrete colours case with no
	     * fermion flow reversal: */

	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		value_type result(0,0);
		Dirac_algebra_type::Id_first(result,it1,it2);
		return result;
	    }
	    
	    /* Contraction algorithm in the continuous colours case with no
	     * fermion flow reversal: */
	    
	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		Dirac_algebra_type::Id_first(result,it1,it2);
	    }

	    /* Contraction algorithm in the discrete colours case with fermion
	     * flow reversal: */

	    static value_type Cc_apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		value_type result(0,0);
		Dirac_algebra_type::Cc_first(result,it1,it2);
		return result;
	    }

	    /* Contraction algorithm in the continuous colours case with fermion
	     * flo reversal: */

	    static void Cc_apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		Dirac_algebra_type::Cc_first(result,it1,it2);
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_SPINOR_CONTR_H_*/

