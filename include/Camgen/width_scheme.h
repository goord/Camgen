//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_WIDTH_SCHEME_H_
#define CAMGEN_WIDTH_SCHEME_H_

#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the width scheme used by the propagators. *
 *                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t>class width_scheme
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);
	
	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);

	public:

	    /// On/off switch for decay widths.

	    static bool switched_on;

	    /// Switch to use complex masses in numerators.

	    static bool complex_masses;
	    
	    /// Switch to use running width scheme.
	    
	    static bool running_widths;

	    /// Switch to use decay widths in spacelike propagators.

	    static bool spacelike_widths;

	    /// Switches decay widths on.

	    static void switch_on()
	    {
		switched_on=true;
	    }

	    /// Switches decay widths off.

	    static void switch_off()
	    {
		switched_on=false;
	    }

	    /// Sets the fixed width scheme for all propagators (no width
	    /// included in spacelike propagators, not gauge-invariant).

	    static void set_fixed_width_scheme()
	    {
		switched_on=true;
		complex_masses=true;
		running_widths=false;
		spacelike_widths=false;
	    }

	    /// Sets the complex mass scheme for all propagators (width included
	    /// in spacelike propagators, allows exact gauge-invariance).
	    
	    static void set_complex_mass_scheme()
	    {
		switched_on=true;
		complex_masses=true;
		running_widths=false;
		spacelike_widths=true;
	    }

	    /// Sets the running-width scheme (no complex masses in numerators,
	    /// no widths in spacelike propagators, not gauge-invariant).

	    static void set_running_width_scheme()
	    {
		switched_on=true;
		complex_masses=false;
		running_widths=true;
		spacelike_widths=false;
	    }

	    /// Evaluates the propagator denominator and complex masses.

	    static void evaluate(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		CAMGEN_ERROR_IF((p==NULL),"NULL momentum pointer encountered");
		
		momentum=p;
		mass=m;
		width=w;
		s=spacetime_type::dot(*p,*p);
		
		if(m==NULL)
		{
		    CAMGEN_MESSAGE_IF((s==0),"On-shell massless propagator encountered");

		    denominator=value_type(0,1)/s;
		    fermion_mass=value_type(0,0);
		    gauge_mass2=s;
		    return;
		}
		if(w==NULL or !switched_on)
		{
		    CAMGEN_MESSAGE_IF((s==(*m)*(*m)),"On-shell zero-width propagator encountered");

		    denominator=value_type(0,1)/(s-(*m)*(*m));
		    fermion_mass=value_type(*m,0);
		    gauge_mass2=value_type((*m)*(*m),0);
		    return;
		}
		if(s<(r_value_type)0 and !spacelike_widths)
		{
		    CAMGEN_MESSAGE_IF((s==(*m)*(*m)),"On-shell spacelike propagator encountered");

		    denominator=value_type(0,1)/(s-(*m)*(*m));
		    fermion_mass=value_type(*m,0);
		    gauge_mass2=value_type((*m)*(*m),0);
		    return;
		}
		if(running_widths)
		{
		    denominator=value_type(0,1)/value_type(s-(*m)*(*m),(*w)*s/(*m));
		}
		else
		{
		    denominator=value_type(0,1)/value_type(s-(*m)*(*m),(*w)*(*m));
		}
		if(complex_masses)
		{
		    gauge_mass2=value_type((*m)*(*m),-(*w)*(*m));
		    fermion_mass=std::sqrt(gauge_mass2);
		}
		else
		{
		    gauge_mass2=value_type((*m)*(*m),0);
		    fermion_mass=value_type(*m,0);
		}
	    }

	protected:

	    /* Momentum flowing through the propagator: */

	    static const momentum_type* momentum;

	    /* Mass of the propagator: */

	    static const r_value_type* mass;

	    /* Width of the propagator: */

	    static const r_value_type* width;

	    /* Invariant mass-squared flowing through the propagator: */

	    static r_value_type s;

	    /* Denominator result from the evaluate() method: */

	    static value_type denominator;

	    /* Complex fermion mass result from the evaluate() method: */

	    static value_type fermion_mass;

	    /* Complex gauge boson mass result from the evaluate() method: */

	    static value_type gauge_mass2;
    };
    template<class model_t>const typename width_scheme<model_t>::momentum_type* width_scheme<model_t>::momentum(NULL);
    template<class model_t>const typename width_scheme<model_t>::r_value_type* width_scheme<model_t>::mass(NULL);
    template<class model_t>const typename width_scheme<model_t>::r_value_type* width_scheme<model_t>::width(NULL);
    template<class model_t>typename width_scheme<model_t>::r_value_type width_scheme<model_t>::s(0);
    template<class model_t>typename width_scheme<model_t>::value_type width_scheme<model_t>::denominator(0,0);
    template<class model_t>typename width_scheme<model_t>::value_type width_scheme<model_t>::fermion_mass(0,0);
    template<class model_t>typename width_scheme<model_t>::value_type width_scheme<model_t>::gauge_mass2(0,0);
    template<class model_t>bool width_scheme<model_t>::switched_on=true;
    template<class model_t>bool width_scheme<model_t>::complex_masses=true;
    template<class model_t>bool width_scheme<model_t>::running_widths=false;
    template<class model_t>bool width_scheme<model_t>::spacelike_widths=true;
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_WIDTH_SCHEME_H_*/

