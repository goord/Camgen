//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_base.h
    \brief Abstract base class for events in Camgen
 */

#ifndef EVT_BASE_H_
#define EVT_BASE_H_

#include <Camgen/sub_proc.h>
#include <Camgen/type_holders.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out>class event_base
    {
        public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef vector<typename model_t::value_type,model_t::dimension> momentum_type;
	    typedef typename momentum_type::value_type value_type;
	    typedef typename momentum_type::size_type size_type;
	    typedef typename get_spacetime_type<model_t>::type spacetime_type;
            typedef particle<model_t> particle_type;
	    
	    static const size_type k0=spacetime_type::timelike_direction;

	    /* Static utility methods: */
	    /*-------------------------*/

	    /// Contracts the argument Lorentz vectors.

	    static value_type dot(const momentum_type& v1,const momentum_type& v2)
	    {
		return spacetime_type::dot(v1,v2);
	    }

	    /// Takes the (Euclidean) inner product of the argument vectors.

	    static value_type space_dot(const momentum_type& v1,const momentum_type& v2)
	    {
		return spacetime_type::space_dot(v1,v2);
	    }

	    /// Takes the invariant mass-squared of the argument vector.

	    static value_type s(const momentum_type& v)
	    {
		return dot(v,v);
	    }

	    /// Takes the invariant mass-squared of the sum of the argument vectors.

	    static value_type s(const momentum_type& v1,const momentum_type& v2)
	    {
		return s(v1+v2);
	    }

	    /// Takes the (signed) invariant mass of the argument vector.

	    static value_type m(const momentum_type& v)
	    {
		value_type sq(s(v));
		return (sq<(value_type)0)?(-std::sqrt(-sq)):(std::sqrt(sq));
	    }

	    /// Takes the (signed) invariant mass of the sum of the argument vectors.

	    static value_type m(const momentum_type& v1,const momentum_type& v2)
	    {
		return m(v1+v2);
	    }

	    /// Takes the energy component of the argument vector.

	    static const value_type& E(const momentum_type& v)
	    {
		return v[k0];
	    }

	    /// Projects out the energy component of the argument vector.

	    static momentum_type vec(const momentum_type& v)
	    {
		momentum_type result(v);
		result[0]=(value_type)0;
		return v;
	    }

	    /// Returns the Euclidean momentum-squared of the argument vector.

	    static value_type P2(const momentum_type& v)
	    {
		return space_dot(v,v);
	    }

	    /// Returns the Euclidean momentum of the argument vector.

	    static value_type P(const momentum_type& v)
	    {
		return std::sqrt(space_dot(v,v));
	    }

	    /// Returns the cosine of the angle between the argument momenta.

	    static value_type cos_alpha(const momentum_type& v1,const momentum_type& v2)
	    {
		return space_dot(v1,v2)/(P(v1)*P(v2));
	    }

	    /// Returns the angle between the argument momenta.

	    static value_type alpha(const momentum_type& v1,const momentum_type& v2)
	    {
		return std::acos(space_dot(v1,v2)/(P(v1)*P(v2)));
	    }

	    /* Default constructor: */

	    event_base():sub_proc(NULL){}

	    /* Destructor: */

	    virtual ~event_base(){}

        private:

            sub_process<model_type,N_in,N_out>* sub_proc;

    };
}

#endif /*EVT_BASE_H_*/
