//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file Lorentz.h
    \brief Lorentz boost functions for Camgen vectors.
 */

#ifndef CAMGEN_LORENTZ_H_
#define CAMGEN_LORENTZ_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Routines involving Lorentz boosts of CAMGEN vectors in D-dimensional *
 * Minkowski spacetime.                                                  *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/vector.h>
#include <Camgen/Minkowski.h>

namespace Camgen
{
    /// Boosts the first vector q, assumed to be defined in the rest-frame of p
    /// and having invariant mass sqrts, to the lab frame.
    
    template<class value_t,std::size_t D>vector<value_t,D>& boost_from_restframe(vector<value_t,D>& q,const vector<value_t,D>& p,const value_t& sqrts)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;

	CAMGEN_ERROR_IF((sqrts!=sqrts),"attempt to boost vector from an undefined frame");
	CAMGEN_ERROR_IF((sqrts==(value_t)0),"attempt to boost vector from a light-like frame");

	value_t pdotq=spacetime_type::space_dot(p,q);
	value_t y=(q[0]+pdotq/(p[0]+sqrts))/sqrts;
	
	q[0]=(p[0]*q[0]+pdotq)/sqrts;
	for(std::size_t mu=1;mu<D;++mu)
	{
	    q[mu]+=y*p[mu];
	}
	return q;
    }
    
    /// Boosts the first vector q, assumed to be defined in the rest-frame of p,
    /// to the lab frame.
    
    template<class value_t,std::size_t D>vector<value_t,D>& boost_from_restframe(vector<value_t,D>& q,const vector<value_t,D>& p)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;

	value_t s=spacetime_type::dot(p,p);
	
	CAMGEN_ERROR_IF((s<(value_t)0),"attempt to boost vector from a space-like frame");
	CAMGEN_ERROR_IF((s==(value_t)0),"attempt to boost vector from a light-like frame");
	
	return boost_from_restframe(q,p,std::sqrt(s));
    }
    
    /// Boosts a copy of the first vector q, assumed to be defined in the rest-frame of p
    /// , assumed to have invariant mass sqrts, to the lab frame.
    
    template<class value_t,std::size_t D>vector<value_t,D> copy_boost_from_restframe(const vector<value_t,D>& q,const vector<value_t,D>& p,const value_t& sqrts)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;

	CAMGEN_ERROR_IF((sqrts!=sqrts),"attempt to boost vector from an undefined frame");
	CAMGEN_ERROR_IF((sqrts==(value_t)0),"attempt to boost vector from a light-like frame");
	
	value_t pdotq=spacetime_type::space_dot(p,q);
	value_t y=(q[0]+pdotq/(p[0]+sqrts))/sqrts;

	vector<value_t,D>qprime;
	qprime[0]=(p[0]*q[0]+pdotq)/sqrts;
	for(std::size_t mu=1;mu<D;++mu)
	{
	    qprime[mu]=q[mu]+y*p[mu];
	}
	return qprime;
    }
    
    /// Boosts a copy of the first vector q, assumed to be defined in the rest-frame of p,
    /// to the lab frame.
    
    template<class value_t,std::size_t D>vector<value_t,D> copy_boost_from_restframe(const vector<value_t,D>& q,const vector<value_t,D>& p)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;

	value_t s=spacetime_type::dot(p,p);
	
	CAMGEN_ERROR_IF((s<(value_t)0),"attempt to boost vector from a space-like frame");
	CAMGEN_ERROR_IF((s==(value_t)0),"attempt to boost vector from a light-like frame");

	return copy_boost_from_restframe(q,p,std::sqrt(s));
    }

    /// Boosts the first vector q to the
    /// rest-frame of the vector p, assumed to have invariant mass sqrts.

    template<class value_t,std::size_t D>vector<value_t,D>& boost_to_restframe(vector<value_t,D>& q,const vector<value_t,D>& p,const value_t& sqrts)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;
	
	CAMGEN_ERROR_IF((sqrts!=sqrts),"attempt to boost vector to an undefined frame");
	CAMGEN_ERROR_IF((sqrts==(value_t)0),"attempt to boost vector to a light-like frame");

	value_t pdotq=spacetime_type::space_dot(p,q);
	value_t y=(q[0]-pdotq/(p[0]+sqrts))/sqrts;
	
	q[0]=(p[0]*q[0]-pdotq)/sqrts;
	for(std::size_t mu=1;mu<D;++mu)
	{
	    q[mu]-=y*p[mu];
	}
	return q;
    }

    /// Boosts the first vector q to the rest-frame of the vector p.
    
    template<class value_t,std::size_t D>vector<value_t,D>& boost_to_restframe(vector<value_t,D>& q,const vector<value_t,D>& p)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;

	value_t s=spacetime_type::dot(p,p);
	
	CAMGEN_ERROR_IF((s<(value_t)0),"attempt to boost vector to a space-like frame");
	CAMGEN_ERROR_IF((s==(value_t)0),"attempt to boost vector to a light-like frame");

	return boost_to_restframe(q,p,std::sqrt(s));
    }
    
    /// Boosts a copy of the first vector q to the
    /// rest-frame of the vector p, assumed to have invariant mass sqrts.

    template<class value_t,std::size_t D>vector<value_t,D> copy_boost_to_restframe(const vector<value_t,D>& q,const vector<value_t,D>& p,const value_t& sqrts)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;
	
	CAMGEN_ERROR_IF((sqrts!=sqrts),"attempt to boost vector to an undefined frame");
	CAMGEN_ERROR_IF((sqrts==(value_t)0),"attempt to boost vector to a light-like frame");
	
	value_t pdotq=spacetime_type::space_dot(p,q);
	value_t y=(q[0]-pdotq/(p[0]+sqrts))/sqrts;

	vector<value_t,D>qprime;
	qprime[0]=(p[0]*q[0]-pdotq)/sqrts;
	for(std::size_t mu=1;mu<D;++mu)
	{
	    qprime[mu]=q[mu]-y*p[mu];
	}
	return qprime;
    }

    /// Boosts a copy of the first vector q to the rest-frame of the vector p.
    
    template<class value_t,std::size_t D>vector<value_t,D> copy_boost_to_restframe(const vector<value_t,D>& q,const vector<value_t,D>& p)
    {
	typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;

	value_t s=spacetime_type::dot(p,p);
	
	CAMGEN_ERROR_IF((s<(value_t)0),"attempt to boost vector to a space-like frame");
	CAMGEN_ERROR_IF((s==(value_t)0),"attempt to boost vector to a light-like frame");

	return copy_boost_to_restframe(q,p,std::sqrt(s));
    }
}

#endif /*CAMGEN_LORENTZ_H_*/


