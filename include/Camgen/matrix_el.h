//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file matrix_el.h
  \brief Abstract base class for user-defined matrix elements.
 */

#ifndef CAMGEN_MATRIX_EL_H_
#define CAMGEN_MATRIX_EL_H_

#include <complex>
#include <Camgen/MC_gen.h>
#include <Camgen/ps_gen_base.h>

namespace Camgen
{
    /// Base class for user-defined matrix elements.

    template<class model_t>class matrix_element: public MC_generator<typename model_t::value_type>,
						 public ps_generator_viewer<model_t>
    {
	public:

	    /* Type definitions: */
	    
	    typedef model_t model_type;
	    typedef ps_generator_viewer<model_t> base_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef std::complex<value_type> c_value_type;
	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::spacetime_type spacetime_type;

	    /// Specifies whether helicity weights must be included.

	    const bool helicity_summed;

	    /// Specifies whether colour weights must be included.

	    const bool colour_summed;

	    /// Constructor, helicity and colour weight included.

	    matrix_element():helicity_summed(false),colour_summed(false){}

	    /// Constructor, where q1 determines whether the given ME is
	    /// helicity-summed and q2 whether it is colour-summed.

	    matrix_element(bool q1,bool q2):helicity_summed(q1),colour_summed(q2){}

	    /// Matrix element evaluation method. Must be implemented by derived
	    /// types.
	    
	    virtual std::complex<value_type> evaluate()=0;

	    /// Generation method. Does nothing, event generation must be
	    /// performed by the external generator.
	    
	    bool generate()
	    {
		return true;
	    }

	    /// Weight evaluation method. Reads in the event from the external
	    /// generator, calls evaluate() to obtain the matrix elements and
	    /// includes phase space, helicity and colour weights, flux factors
	    /// etc.

	    bool evaluate_weight()
	    {
		if(this->gen!=NULL)
		{
		    this->weight()=this->gen->ps_weights(!helicity_summed,!colour_summed);
		    if(this->weight()>(value_type)0)
		    {
			this->integrand()=this->gen->ps_factors(true,true)*std::norm(evaluate());
		    }
		    return true;
		}
		else
		{
		    CAMGEN_MESSAGE("Phase space generator was not defined in matrix element class...setting the weight zero");
		    return false;
		}
	    }

	    /// Returns the phase space integration prefactors, obtained from
	    /// the external generator.

	    value_type ps_factors()
	    {
		if(this->gen==NULL)
		{
		    return 0;
		}
		return this->gen->ps_factors(true,true);
	    }

	    /// Returns the full weight.

	    value_type w() const
	    {
		return this->weight()*this->integrand();
	    }

	private:

	    bool auto_update;
    };

    /// Unit matrix element, useful for computing the phase space volume.

    template<class model_t>class unit_matrix_element: public matrix_element<model_t>
    {
	public:

	    /* Type definitions: */

	    typedef matrix_element<model_t> base_type;
	    typedef typename model_t::value_type value_type;
	    
	    /// Constructor.

	    unit_matrix_element(){}

	    /// Matrix element evaluation implementation. Returns unity.
	    
	    virtual std::complex<value_type> evaluate()
	    {
		return std::complex<value_type>(1,0);
	    }
    };

    /// Reweighting matrix element.

    template<class model_t,class model_t2>class matrix_element_reweighter: public matrix_element<model_t>
    {
	public:

	    /* Type definitions: */

	    typedef matrix_element<model_t> base_type;
	    typedef typename model_t::value_type value_type;
	    
	    /// Constructor.
	    //
	    // 
	    //
	    //
    };
}

#endif /*CAMGEN_MATRIX_EL_H_*/

