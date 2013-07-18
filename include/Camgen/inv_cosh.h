//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_INV_COSH_H_
#define CAMGEN_INV_COSH_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the rapidity generator                *
 *                                                                     *
 *                          g(y)~1/cosh(y-y0)                          *
 *                                                                     *
 * where y0 denotes the central maximum of the generated distribution. *
 *                                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <sstream>
#include <Camgen/plt_strm.h>
#include <Camgen/inv_gen.h>

namespace Camgen
{
    template<class value_t,class rng_t>class inv_cosh_y_generator: public inversion_s_generator<value_t,rng_t,inv_cosh_y_generator<value_t,rng_t> >
    {
	typedef inversion_s_generator< value_t,rng_t,inv_cosh_y_generator<value_t,rng_t> > direct_base_type;
	typedef s_generator<value_t,rng_t> base_type;

	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Public data: */

	    const value_type* const y0;

	    /* Default constructor: */

	    inv_cosh_y_generator(const value_type* y0_=NULL):y0(y0_)
	    {
		refresh_params();
		this->refresh_s_min();
		this->refresh_s_max();
	    }

	    /* Constructor with central value argument: */

	    inv_cosh_y_generator(value_type* y_,const value_type* y0_=NULL):direct_base_type(y_),y0(y0_)
	    {
		refresh_params();
		this->refresh_s_min();
		this->refresh_s_max();
	    }

	    /* Copy constructor: */

	    inv_cosh_y_generator(const inv_cosh_y_generator<value_t,rng_t>& other):direct_base_type(other),y0(other.y0){}

	    /* Public modifiers: */
	    /*-------------------*/

	    bool refresh_params()
	    {
		this->min_cumulant_lim=(value_type)0;
		this->max_cumulant_lim=pi;
		return (this->refresh_s_min() and this->refresh_s_max());

	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method: */

	    inv_cosh_y_generator<value_t,rng_t>* clone() const
	    {
		return new inv_cosh_y_generator<value_t,rng_t>(*this);
	    }

	    /* Central peak readout: */

	    value_type peak() const
	    {
		return (y0==NULL)?(value_type)0:(*y0);
	    }

	    /* Probability function: */

	    value_type pdf(const value_type& y) const
	    {
		value_type x(std::exp(y-peak()));
		return (value_type)2*x/(x*x+(value_type)1);
	    }

	    /* Cumulant: */

	    value_type cdf(const value_type& y) const
	    {
		return (value_type)2*std::atan(std::exp(y-peak()));
	    }

	    /* Inverse of the cumulant: */

	    value_type inverse_cdf(const value_type rho) const
	    {
		return (peak()+std::log(std::tan(0.5*rho)));
	    }
	    
	    /* Double-dispatch integrator functions. */
	    
	    value_type integrate_with(const s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		log(log_level::error)<<CAMGEN_STREAMLOC<<"integration method not defined"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const BW_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		log(log_level::error)<<CAMGEN_STREAMLOC<<"integration method not defined"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const pl_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		log(log_level::error)<<CAMGEN_STREAMLOC<<"integration method not defined"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const uni_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		log(log_level::error)<<CAMGEN_STREAMLOC<<"integration method not defined"<<endlog;
		return 0;
	    }
	    value_type integrate_with(const Dd_s_generator<value_t,rng_t>* gen,const value_type& sqrts) const
	    {
		log(log_level::error)<<CAMGEN_STREAMLOC<<"integration method not defined"<<endlog;
		return 0;
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Overridden printing function: */

	    std::ostream& print(std::ostream& os) const
	    {
		os<<"1/cosh("<<peak()<<')';
		return os;
	    }

	    /* Outputs the type of the generator. */

	    std::string type() const
	    {
		return "invcosh";
	    }

	    /* Outputs the detailed type of the generator. */

	    std::string detailed_type() const
	    {
		std::stringstream ss;
		ss<<"1/cosh("<<peak()<<')';
		return ss.str();
	    }

	    /* Plotting utility method. */

	    function_stream* pdf_plot(const value_type& x) const
	    {
		value_type n=x/this->norm;
		std::ostringstream sstrm;
		sstrm<<n<<"/cosh(x-"<<peak()<<')';
		function_stream* result=new function_stream(sstrm.str());
		result->title="density";
		result->style="lines";
		return result;
	    }

	private:

	    static const value_t pi;
    };
    template<class value_t,class rng_t>const value_t inv_cosh_y_generator<value_t,rng_t>::pi(std::acos(-(value_t)1));
}

#endif /*CAMGEN_INV_COSH_H_*/

