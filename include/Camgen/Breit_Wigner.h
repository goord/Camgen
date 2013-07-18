//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_BREIT_WIGNER_H_
#define CAMGEN_BREIT_WIGNER_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the Breit-Wigner invariant mass generator *
 *                                                                         *
 *                         g(s)~1/((s-m^2)^2+m^2*w^2)                      *
 *                                                                         *
 * where m denotes the mass and w the width (fixed) of the resonance.      *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <complex>
#include <sstream>
#include <Camgen/plt_strm.h>
#include <Camgen/inv_gen.h>

namespace Camgen
{
    /* Momentum channel class template generating invariant masses according to
     * a Breit-Wigner distribution. */

    template<class value_t,class rng_t>class BW_s_generator: public inversion_s_generator< value_t,rng_t,BW_s_generator<value_t,rng_t> >
    {
	friend class s_integrator<value_t,rng_t>;
	typedef inversion_s_generator< value_t,rng_t,BW_s_generator<value_t,rng_t> > direct_base_type;
	typedef s_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Public data: */
	    /*--------------*/

	    /* Mass and width addresses: */

	    const value_type* const m;
	    const value_type* const w;

	    /* Public constuctors: */
	    /*---------------------*/
	    
	    /* Default constructor: */

	    BW_s_generator(const value_type* m_,const value_type* w_):m(m_),w(w_)
	    {
		refresh_params();
		this->refresh_s_min();
		this->refresh_s_max();
	    }
	    
	    /* Constructor with invariant mass pointer argument. */
	    
	    BW_s_generator(value_type* s_,const value_type* m_,const value_type* w_):direct_base_type(s_),m(m_),w(w_)
	    {
		refresh_params();
		this->refresh_s_min();
		this->refresh_s_max();
	    }

	    /* Copy constructor: */

	    BW_s_generator(const BW_s_generator<value_t,rng_t>& other):direct_base_type(other),m(other.m),w(other.w),mm(other.mm),mw(other.mw),mw2(other.mw2),alpha(other.alpha),xi_plus(other.xi_plus),xi_min(other.xi_min),nu(other.nu),nubar(other.nubar),lognu(other.lognu),logminnu(other.logminnu),lognubar(other.lognubar),logminnubar(other.logminnubar){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Updates the internal data from the current values of the mass
	    and width pointers. */

	    bool refresh_params()
	    {
		mm=(m==NULL)?(value_type)0:((*m)*(*m));
		mw=(m==NULL)?(value_type)0:((w==NULL)?(value_type)0:((*m)*(*w)));
		if(mw==(value_type)0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid Breit-Wigner parameters mass: "<<mass()<<" GeV, width: "<<width()<<" GeV encountered"<<endlog;
		    return false;
		}
		this->max_cumulant_lim=pi2/mw;
		this->min_cumulant_lim=-this->max_cumulant_lim;
		mw2=mw*mw;
		alpha=std::sqrt(mm*mm+mw2);
		xi_plus=std::sqrt((value_type)2*(alpha+mm));
		xi_min=std::sqrt((value_type)2*(alpha-mm));
		nu=std::sqrt(std::complex<value_type>(mm,-mw));
		nubar=std::sqrt(std::complex<value_type>(mm,mw));
		return (this->refresh_s_min() and this->refresh_s_max());
	    }
	    
	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method: */

	    BW_s_generator<value_t,rng_t>* clone() const
	    {
		return new BW_s_generator<value_t,rng_t>(*this);
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
	
	    /* Returns the width value. */

	    value_type width() const
	    {
		return (w==NULL)?(value_type)0:(*w);
	    }

	    /* Probability density function. */

	    value_type pdf(const value_type& s) const
	    {
		return (value_type)1/((s-mm)*(s-mm)+mw2);
	    }

	    /* Cumulator distribution function. */

	    value_type cdf(const value_type& s) const
	    {
		return std::atan((s-mm)/mw)/mw;
	    }

	    /* Inverse of the cumulator. */

	    value_type inverse_cdf(const value_type& rho) const
	    {
		return (mm+mw*std::tan(mw*rho));
	    }
	    
	    /* Double-dispatch integrator functions. */
	    
	    value_type integrate_with(const s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return gen->integrate_with(this,sqrts);
	    }
	    value_type integrate_with(const BW_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(this,gen,sqrts);
	    }
	    value_type integrate_with(const pl_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		typedef pl_s_generator<value_t,rng_t> gen_t;
		return s_integrator<value_t,rng_t>::template integrate<gen_t>(this,static_cast<const inversion_s_generator<value_t,rng_t,gen_t>*>(gen),sqrts);
	    }
	    value_type integrate_with(const uni_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		return s_integrator<value_t,rng_t>::integrate(this,gen,sqrts);
	    }
	    value_type integrate_with(const Dd_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		typedef BW_s_generator<value_t,rng_t> gen_t;
		return s_integrator<value_t,rng_t>::template integrate<gen_t>(static_cast<const direct_base_type*>(this),gen,sqrts);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Overridden printing function: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<"BW("<<mass()<<','<<width()<<')';
		return os;
	    }

	    /* Outputs the type of the generator. */

	    std::string type() const
	    {
		return "BW";
	    }

	    /* Outputs the detailed type of the generator. */

	    std::string detailed_type() const
	    {
		std::stringstream ss;
		ss<<"BW("<<mass()<<','<<width()<<')';
		return ss.str();
	    }

	    /* Plotting utility method. */

	    function_stream* pdf_plot(const value_type& x) const
	    {
		value_type n=x/this->norm;
		std::ostringstream sstrm;
		sstrm<<n<<"/((x-"<<mm<<")**2+"<<mw2<<")";
		function_stream* result=new function_stream(sstrm.str());
		result->title="density";
		result->style="lines";
		return result;
	    }

	private:

	    /* Extra internal parameters: */

	    value_type mm,mw,mw2,alpha,xi_plus,xi_min;
	    std::complex<value_type>nu,nubar,lognu,logminnu,lognubar,logminnubar;

	    /* One-half pi: */

	    static const value_t pi2;
    };
    template<class value_t,class rng_t>const value_t BW_s_generator<value_t,rng_t>::pi2(std::acos((value_t)0));
}

#endif /*CAMGEN_BREIT_WIGNER_H_*/

