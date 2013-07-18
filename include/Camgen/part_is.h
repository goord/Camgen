//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file part_is.h
    \brief Implementation of fixed CM-energy 1- and 2-particle partonic initial state generators.
 */

#ifndef CAMGEN_PART_IS_H_
#define CAMGEN_PART_IS_H_

#include <Camgen/Minkowski.h>
#include <Camgen/init_state.h>
#include <Camgen/MC_config.h>
#include <Camgen/ps_vol.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Standard initial state generator classes. Implementation are only provided  *
 * for 1- and 2-particle initial states, and the momenta are always created in *
 * the CM frame, with an effective invariant mass equal to the total invariant *
 * mass.                                                                       *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Standard 1-particle initial-state implementation: */

    template<class model_t>class partonic_is<model_t,1,Minkowski_type>: public initial_state<model_t,1>
    {
	/* Type definitions: */

	typedef initial_state<model_t,1> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::size_type size_type;

	    /* Static integer constants: */
	    
	    static const size_type k0=spacetime_type::timelike_direction;
	    static const size_type kL=(model_type::beam_direction<0)?(-model_type::beam_direction):(model_type::beam_direction);
	    static const int kL_sign=(model_type::beam_direction<0)?(-1):1;

	    /* Public constructors */
	    /*---------------------*/

	    /* Default constructor. */
	    
	    partonic_is():base_type(false,false,true){}

	    /* Constructor with beam energy. */
	    
	    partonic_is(const value_type* E_beam):base_type(base_type::make_vector(E_beam),false,false,true){}

	    /* Copy constructor. */

	    partonic_is(const partonic_is<model_t,1,Minkowski_type>& other):base_type(other){}

	    /* Public modifiers */
	    /*------------------*/

	    /* Hadronic invariant mass refresher. */

	    bool refresh_Ecm()
	    {
		return (this->set_Ecm(this->m(0)) and this->set_Ecm_hat(this->m(0)));
	    }

	    /* Generation implementation. Ignores any invariant mass generation
	    in the output channel. */

	    bool generate()
	    {
		if(!this->set_Ecm_hat(this->Ecm()))
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		momentum_type& p=this->p(0);
		p.assign(0);
		value_type M=this->m(0);
		value_type E=this->beam_energy(0);
		if(E==(value_type)0)
		{
		    p[k0]=M;
		    this->weight()=(value_type)1;
		    return true;
		}
		p[k0]=E+M;
		p[kL]=kL_sign*std::sqrt(E*(E+(value_type)2*M));
		this->weight()=(value_type)1;
		return true;
	    }

	    bool generate_p()
	    {
		momentum_type& p=this->p(0);
		p.assign(0);
		value_type M=this->m(0);
		value_type E=this->beam_energy(0);
		if(E==(value_type)0)
		{
		    p[k0]=M;
		    this->weight()=(value_type)1;
		    return true;
		}
		p[k0]=E+M;
		p[kL]=k0*std::sqrt(E*(E+(value_type)2*M));
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Returns the flux factor. */
	    
	    value_type flux_factor() const
	    {
		return (value_type)0.5/this->m(0);
	    }

	    /* Weight evaluation method. */

	    bool evaluate_weight()
	    {
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /* Clone method implementation. */

	    partonic_is<model_t,1,Minkowski_type>* clone() const
	    {
		return new partonic_is<model_t,1,Minkowski_type>(*this);
	    }

	    /* Beam particle id: */

	    int beam_id(size_type i) const
	    {
		return (i==0)?(this->parton_id(i)):0;
	    }

	    /* cernlib pdf group number: */

	    int pdfg(size_type i) const
	    {
		return (i==0)?(-1):0;
	    }

	    /* Returns the cernlib pdf set number: */

	    int pdfs(size_type i) const
	    {
		return (i==0)?(-1):0;
	    }

	    /* Serialisation: */
	    /*----------------*/

	    std::string type() const
	    {
		return "partonic";
	    }
    };
    template<class model_t>const typename partonic_is<model_t,1,Minkowski_type>::size_type partonic_is<model_t,1,Minkowski_type>::k0;
    template<class model_t>const typename partonic_is<model_t,1,Minkowski_type>::size_type partonic_is<model_t,1,Minkowski_type>::kL;
    template<class model_t>const int partonic_is<model_t,1,Minkowski_type>::kL_sign;

    /* Standard 2-particle initial-state implementation. */

    template<class model_t>class partonic_is<model_t,2,Minkowski_type>: public initial_state<model_t,2>
    {
	/* Type definition: */

	typedef initial_state<model_t,2> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::size_type size_type;

	    /* Static integer constants: */
	    
	    static const size_type k0=spacetime_type::timelike_direction;
	    static const size_type kL=(model_type::beam_direction<0)?(-model_type::beam_direction):(model_type::beam_direction);
	    static const int kL_sign=(model_type::beam_direction<0)?(-1):1;

	    /* Public constructors */
	    /*---------------------*/

	    /* Constructor with no beam energies. */
	    
	    partonic_is():base_type(false,false,true){}

	    /* Constructor with equal beam energies. */
	    
	    partonic_is(const value_type& E_beam):base_type(false,false,true)
	    {
		this->E_beams[0]=(value_type)0.5*E_beam;
		this->E_beams[1]=(value_type)0.5*E_beam;
	    }

	    /* Constructor with distinguished beam energies. */
	    
	    partonic_is(const value_type& E_beam1,const value_type& E_beam2):base_type(false,false,true)
	    {
		this->E_beams[0]=E_beam1;
		this->E_beams[1]=E_beam2;
	    }

	    /* Copy constructor. */

	    partonic_is(const partonic_is<model_t,2,Minkowski_type>& other):base_type(other){}

	    /* Public modifiers */
	    /*------------------*/

	    /* Hadronic invariant mass refresher. */

	    bool refresh_Ecm()
	    {
		value_type M1=this->m(0);
		value_type M2=this->m(1);
		value_type E1=this->beam_energy(0)+M1;
		value_type E2=this->beam_energy(1)+M2;
		value_type P1=std::sqrt(E1*E1-M1*M1);
		value_type P2=std::sqrt(E2*E2-M2*M2);
		return (this->set_s(M1*M1+M2*M2+(value_type)2*(E1*E2+P1*P2)) and this->set_Ecm_hat(this->Ecm()));
	    }

	    /* Generation implementation. */

	    bool generate()
	    {
		if(!this->set_Ecm_hat(this->Ecm()))
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->p(0).assign(0);
		this->p(1).assign(0);
		value_type M1=this->m(0);
		value_type M2=this->m(1);
		value_type E1=this->beam_energy(0)+M1;
		value_type E2=this->beam_energy(1)+M2;
		value_type P1=std::sqrt(E1*E1-M1*M1);
		value_type P2=std::sqrt(E2*E2-M2*M2);
		this->p(0)[k0]=E1;
		this->p(0)[kL]=kL_sign*P1;
		this->p(1)[k0]=E2;
		this->p(1)[kL]=-kL_sign*P2;
		this->weight()=(value_type)1;
		return true;
	    }

	    bool generate_p()
	    {
		this->p(0).assign(0);
		this->p(1).assign(0);
		value_type M1=this->m(0);
		value_type M2=this->m(1);
		value_type E1=this->beam_energy(0)+M1;
		value_type E2=this->beam_energy(1)+M2;
		value_type P1=std::sqrt(E1*E1-M1*M1);
		value_type P2=std::sqrt(E2*E2-M2*M2);
		this->p(0)[k0]=E1;
		this->p(0)[kL]=kL_sign*P1;
		this->p(1)[k0]=E2;
		this->p(1)[kL]=-kL_sign*P2;
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Returns the flux factor. */
	    
	    value_type flux_factor() const
	    {
		value_type E1=this->p(0)[k0];
		value_type E2=this->p(1)[k0];
		value_type P1=std::abs(this->p(0)[kL]);
		value_type P2=std::abs(this->p(1)[kL]);
		return (value_type)0.25/(E1*P2+E2*P1);
	    }

	    /* Weight evaluation method. */

	    bool evaluate_weight()
	    {
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /* Clone method implementation. */

	    partonic_is<model_t,2,Minkowski_type>* clone() const
	    {
		return new partonic_is<model_t,2,Minkowski_type>(*this);
	    }

	    /* Beam particle id: */

	    int beam_id(size_type i) const
	    {
		if(i<2)
		{
		    return this->parton_id(i);
		}
		return 0;
	    }

	    /* cernlib pdf group number: */

	    int pdfg(size_type i) const
	    {
		return (i<2)?(-1):0;
	    }

	    /* Returns the cernlib pdf set number: */

	    int pdfs(size_type i) const
	    {
		return (i<2)?(-1):0;
	    }

	    /* Serialisation: */
	    /*----------------*/

	    std::string type() const
	    {
		return "partonic";
	    }
    };
    template<class model_t>const typename partonic_is<model_t,2,Minkowski_type>::size_type partonic_is<model_t,2,Minkowski_type>::k0;
    template<class model_t>const typename partonic_is<model_t,2,Minkowski_type>::size_type partonic_is<model_t,2,Minkowski_type>::kL;
    template<class model_t>const int partonic_is<model_t,2,Minkowski_type>::kL_sign;
}

#endif /*CAMGEN_PART_IS_H_*/

