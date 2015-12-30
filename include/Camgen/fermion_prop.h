//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FERMION_PROP_H_
#define CAMGEN_FERMION_PROP_H_

#include <vector>
#include <Camgen/width_scheme.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of massless and massive fermion and anti-fermion *
 * propagator classes.                                                         *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Fermion propagator: */

    template<class model_t>class fermion_propagator: public width_scheme<model_t>
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime and Dirac algebra type definitions: */

	    DEFINE_SPACETIME_TYPE(model_t);
	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	public:

	    /* Size of a propagating subamplitude tensor: */

	    static const size_type size=Dirac_algebra_type::base_type::index_range;

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
		Dirac_algebra_type::initialise();
	    }

	    /* Function returning the index ranges of propagating (colourless)
	     * subamplitude tensors: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& v)
	    {
		v.clear();
		v.push_back(Dirac_algebra_type::base_type::index_range);
		return v;
	    }

	    /* Evaluates the propagator denominator: */

	    static void refresh(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		width_scheme<model_t>::evaluate(p,m,w);
	    }

	    /* Function propagating the spinor 'first' points to by calling the
	     * appropriate member function in the Dirac algebra class: */

	    static void evaluate(iterator first)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		
		if(width_scheme<model_t>::mass==NULL)
		{
		    Dirac_algebra_type::massless_prop(width_scheme<model_t>::denominator,first,*p);
		    return;
		}
		Dirac_algebra_type::massive_prop(width_scheme<model_t>::denominator,width_scheme<model_t>::fermion_mass,first,*p);
	    }

	    /* Function propagating all spinors between 'first' and 'last' by
	     * calling the appropriate member function in the Dirac algebra
	     * class: */

	    static void evaluate_range(iterator first,iterator last)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		
		if(width_scheme<model_t>::mass==NULL)
		{
		    Dirac_algebra_type::massless_prop(width_scheme<model_t>::denominator,first,last,*p);
		    return;
		}
		Dirac_algebra_type::massive_prop(width_scheme<model_t>::denominator,width_scheme<model_t>::fermion_mass,first,last,*p);
	    }
    };
    template<class model_t>const typename fermion_propagator<model_t>::size_type fermion_propagator<model_t>::size;

    /* Anti-fermion propagator: */

    template<class model_t>class anti_fermion_propagator: public width_scheme<model_t>
    {
	public:

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime and Dirac algebra type definitions: */

	    DEFINE_SPACETIME_TYPE(model_t);
	    DEFINE_DIRAC_ALGEBRA_TYPE(model_t);

	public:	

	    /* Size of a propagating subamplitude tensor: */

	    static const size_type size=Dirac_algebra_type::base_type::index_range;

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
		Dirac_algebra_type::initialise();
	    }

	    /* Function returning the index ranges of propagating (colourless)
	     * subamplitude tensors: */

	    static std::vector<size_type>& get_index_ranges(std::vector<size_type>& v)
	    {
		v.clear();
		v.push_back(Dirac_algebra_type::base_type::index_range);
		return v;
	    }

	    /* Evaluates the propagator denominator: */

	    static void refresh(const momentum_type* p,const r_value_type* m,const r_value_type* w)
	    {
		width_scheme<model_t>::evaluate(p,m,w);
	    }

	    /* Function propagating the anti-spinor 'first' points to by calling the
	     * appropriate member function in the Dirac algebra class: */

	    static void evaluate(iterator first)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");
		
		if(width_scheme<model_t>::mass==NULL)
		{
		    Dirac_algebra_type::massless_anti_prop(width_scheme<model_t>::denominator,first,*p);
		    return;
		}
		Dirac_algebra_type::massive_anti_prop(width_scheme<model_t>::denominator,width_scheme<model_t>::fermion_mass,first,*p);
	    }

	    /* Function propagating all anti-spinors between 'first' and 'last' by
	     * calling the appropriate member function in the Dirac algebra
	     * class: */

	    static void evaluate_range(iterator first,iterator last)
	    {
		const momentum_type* p=width_scheme<model_t>::momentum;
		CAMGEN_ERROR_IF((p==NULL),"attempt to dereference a NULL momentum pointer");

		if(width_scheme<model_t>::mass==NULL)
		{
		    Dirac_algebra_type::massless_anti_prop(width_scheme<model_t>::denominator,first,last,*p);
		    return;
		}
		Dirac_algebra_type::massive_anti_prop(width_scheme<model_t>::denominator,width_scheme<model_t>::fermion_mass,first,last,*p);
	    }
    };
    template<class model_t>const typename anti_fermion_propagator<model_t>::size_type anti_fermion_propagator<model_t>::size;
}

#include <Camgen/undef_args.h>

#endif /*FERMION_PROP_H_*/

