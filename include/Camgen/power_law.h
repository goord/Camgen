//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_POWER_LAW_H_
#define CAMGEN_POWER_LAW_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the invariant mass generator according to       *
 *                                                                               *
 * 				g(s)~1/|s-m^2|^nu                                *
 *                                                                               *
 * where nu can be any number. If larger then 1 however, the density contains a  *
 * non-integrable pole at s=m^2 and bounds should be chosen not to contain this  *
 * singularity.                                                                  *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <sstream>
#include <Camgen/plt_strm.h>
#include <Camgen/inv_gen.h>

namespace Camgen
{
    /* Momentum channel class template generating invariant masses according to
     * a power-law distribution. */
    
    template<class value_t,class rng_t>class pl_s_generator: public inversion_s_generator< value_t,rng_t,pl_s_generator<value_t,rng_t> >
    {
	friend class s_integrator<value_t,rng_t>;
	typedef inversion_s_generator< value_t,rng_t,pl_s_generator<value_t,rng_t> > direct_base_type;
	typedef s_generator<value_t,rng_t> base_type;

	public:

	    enum exponent_type
	    {
		positive,
		zero,
		negative_integrable,
		logarithmic,
		negative
	    };

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Public data members: */
	    /*----------------------*/

	    /* Offset and exponent adresses: */

	    const value_type* const m;
	    const value_type* const nu;

	    /* Constructors: */
	    /*---------------*/
	    
	    /* Default constructor: */

	    pl_s_generator(const value_type* m_,const value_type* nu_):m(m_),nu(nu_)
	    {
		refresh_params();
		this->refresh_s_min();
		this->refresh_s_max();
	    }
	    
	    /* Constructor with inveriant mass pointer argument. */
	    
	    pl_s_generator(value_type* s_,const value_type* m_,const value_type* nu_):direct_base_type(s_),m(m_),nu(nu_)
	    {
		refresh_params();
		this->refresh_s_min();
		this->refresh_s_max();
	    }

	    /* Copy constructor: */

	    pl_s_generator(const pl_s_generator<value_t,rng_t>& other):direct_base_type(other),m(other.m),nu(other.nu),exp_t(other.exp_t),nu_inv(other.nu_inv),mm(other.mm){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Updates the internal data from the current values of the mass and
	     * width pointers. */

	    bool refresh_params()
	    {
		if(nu==NULL)
		{
		    exp_t=zero;
		    this->min_cumulant_lim=(value_type)0;
		    this->max_cumulant_lim=std::numeric_limits<value_type>::infinity();
		}
		else
		{
		    if(*nu>(value_type)1)
		    {
			exp_t=negative;
			this->min_cumulant_lim=(value_type)0;
			this->max_cumulant_lim=(value_type)0;
		    }
		    if(*nu==(value_type)1)
		    {
			exp_t=logarithmic;
			this->min_cumulant_lim=(value_type)0;
			this->max_cumulant_lim=std::numeric_limits<value_type>::infinity();
		    }
		    if(*nu<(value_type)1 and *nu>(value_type)0)
		    {
			exp_t=negative_integrable;
			this->min_cumulant_lim=(value_type)0;
			this->max_cumulant_lim=std::numeric_limits<value_type>::infinity();
		    }
		    if(*nu==(value_type)0)
		    {
			exp_t=zero;
			this->min_cumulant_lim=(value_type)0;
			this->max_cumulant_lim=std::numeric_limits<value_type>::infinity();
		    }
		    if(*nu<(value_type)0)
		    {
			exp_t=positive;
			this->min_cumulant_lim=(value_type)0;
			this->max_cumulant_lim=std::numeric_limits<value_type>::infinity();
		    }
		}
		nu_inv=(value_type)1/((value_type)1-exp());
		mm=(m==NULL)?(value_type)0:((*m)*(*m));
		return (this->refresh_s_min() and this->refresh_s_max());
	    }

	    /* Public const methods: */
	    /*-----------------------*/
	    
	    /* Clone method implementation: */

	    pl_s_generator<value_t,rng_t>* clone() const
	    {
		return new pl_s_generator<value_t,rng_t>(*this);
	    }

	    /* Returns the mass value. */

	    value_type mass() const
	    {
		return (m==NULL)?(value_type)0:(*m);
	    }

	    /* Returns the mass-squared value. */

	    value_type mass2() const
	    {
		return mm;
	    }
	
	    /* Returns the exponent value. */

	    value_type exp() const
	    {
		return (nu==NULL)?(value_type)0:(*nu);
	    }

	    /* Probability density function. */

	    value_type pdf(const value_type& s) const
	    {
		return (nu==NULL)?(value_type)1:std::pow(std::abs(s-mm),-*nu);
	    }

	    /* Cumulator distribution function. */

	    value_type cdf(const value_type& s) const
	    {
		switch(exp_t)
		{
		    case zero:
			return s;
		    case logarithmic:
			return (s<mm)?(-std::log(mm-s)):std::log(s-mm);
		    default:
			return (s<mm)?(-nu_inv*std::pow(mm-s,(value_type)1-*nu)):(nu_inv*std::pow(s-mm,(value_type)1-*nu));
		}
	    }

	    /* Inverse of the cumulator. */

	    value_type inverse_cdf(const value_type& rho) const
	    {
		switch(exp_t)
		{
		    case zero:
			return rho;
		    case logarithmic:
			return ((this->s_min()>mm)?(mm+std::exp(rho)):((this->s_max()<mm)?(mm-std::exp(-rho)):mm));
		    default:
			return ((rho/nu_inv)>(value_type)0)?(mm+std::pow(rho/nu_inv,nu_inv)):(mm-std::pow(-rho/nu_inv,nu_inv));
		}
	    }

	    /* Double-dispatch integrator functions. */
	    
	    value_type integrate_with(const s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return gen->integrate_with(this,sqrts);
	    }
	    value_type integrate_with(const BW_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(gen,this,sqrts);
	    }
	    value_type integrate_with(const pl_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		typedef pl_s_generator<value_t,rng_t> type;
		return s_integrator<value_t,rng_t>::template integrate<type,type>(static_cast<const direct_base_type*>(this),static_cast<const direct_base_type*>(gen),sqrts);
	    }
	    value_type integrate_with(const uni_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(this,gen,sqrts);
	    }
	    value_type integrate_with(const Dd_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		typedef pl_s_generator<value_t,rng_t> type;
		return s_integrator<value_type,rng_t>::template integrate<type>(static_cast<const direct_base_type*>(this),gen,sqrts);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Printing method: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<"pl("<<mass()<<','<<exp()<<')';
		return os;
	    }

	    /* Outputs the type of the generator. */

	    std::string type() const
	    {
		return "pl";
	    }

	    /* Outputs the detailed type of the generator. */

	    std::string detailed_type() const
	    {
		std::stringstream ss;
		ss<<"pl("<<mass()<<','<<exp()<<')';
		return ss.str();
	    }

	    /* Plotting utility method. */

	    function_stream* pdf_plot(const value_type& x) const
	    {
		value_type n=x/this->norm;
		std::ostringstream sstrm;
		switch(exp_t)
		{
		    case positive:
			sstrm<<n<<"*abs(x-"<<mm<<")**"<<(-*nu);
			break;
		    case zero:
			sstrm<<n;
			break;
		    case negative_integrable:
			sstrm<<n<<"/abs(x-"<<mm<<")**"<<*nu;
			break;
		    case logarithmic:
			sstrm<<n<<"/abs(x-"<<mm<<")";
			break;
		    default:
			sstrm<<n<<"/abs(x-"<<mm<<")**"<<*nu;
		}
		function_stream* result=new function_stream(sstrm.str());
		result->title="density";
		result->style="lines";
		return result;
	    }

	protected:

	    /* Normalisation refresher. */

	    void refresh_norm()
	    {
		if((exp_t==logarithmic or exp_t==negative) and mm>this->s_min() and mm<this->s_max())
		{
		    this->norm=std::numeric_limits<value_type>::infinity();
		}
		else
		{
		    this->norm=std::abs(this->max_cumulant-this->min_cumulant);
		}
	    }

	private:

	    /* Extra parameters: */

	    exponent_type exp_t;
	    value_type nu_inv,mm;
    };
}

#endif /*CAMGEN_POWER_LAW_H_*/

