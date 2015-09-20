//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file val_gen_fac.h
    \brief Factory for value Monte Carlo generators.
 */

#ifndef CAMGEN_VAL_GEN_FAC_H_
#define CAMGEN_VAL_GEN_FAC_H_

#include <Camgen/val_gen.h>
#include <Camgen/val_gen_grid.h>
#include <Camgen/inv_cosh.h>
#include <Camgen/Dirac_delta.h>
#include <Camgen/Breit_Wigner.h>
#include <Camgen/power_law.h>
#include <Camgen/uni_val_gen.h>

/* * * * * * * * * * * * * * * * * * * * * * * *
 * Factory for invariant mass sampler classes. *
 *                                             *
 * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class value_t,class rng_t>class value_generator_factory
    {
	public:

	    /* Type definitions: */
	    
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    typedef value_generator<value_t,rng_t> value_generator_type;
	    typedef std::size_t size_type;

	    /// Creates an invariant mass sampler depending on given parameters.

	    static value_generator_type* create_instance(const value_type* mass, const value_type* width, const value_type* nu,bool adaptive=false)
	    {
		value_generator_type* inversion_generator=NULL;
		bool on_shell=false;

		if(width==NULL and nu==NULL)
		{
		    inversion_generator=new Dirac_delta<value_type,rng_t>(mass);
		    on_shell=true;
		}
		if(width!=NULL)
		{
		    inversion_generator=new Breit_Wigner<value_type,rng_t>(mass,width);
		}
		else if(nu!=NULL)
		{
		    inversion_generator=new power_law<value_type,rng_t>(mass,nu);
		}

		if(adaptive and !on_shell)
		{
		    return new adaptive_value_generator<value_type,rng_t>(inversion_generator,grid_bins(),grid_mode());
		}

		return inversion_generator;
	    }
    };
}

#endif /*CAMGEN_S_GEN_FAC_H_*/
