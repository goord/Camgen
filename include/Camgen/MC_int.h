//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file MC_int.h
    \brief Abstract base class template for Monte Carlo integrator classes.
 */

#ifndef CAMGEN_MC_INT_H_
#define CAMGEN_MC_INT_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for the Monte Carlo integrators *
 *                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/num_utils.h>
#include <Camgen/logstream.h>
#include <Camgen/MC_int.h>
#include <Camgen/whist.h>
#include <Camgen/rn_strm.h>
#include <Camgen/MC_integral.h>
#include <Camgen/MC_gen.h>
#include <Camgen/MC_int_base.h>

namespace Camgen
{
    template<class value_t>class MC_integrator: public MC_generator<value_t>, public MC_integrator_base<value_t>
    {
	public:

	    /* Public type definitions: */
	    /*--------------------------*/

	    typedef std::size_t size_type;

	    /// Value type definition.
	    
	    typedef value_t value_type;

	    /// Evaluated Monte Carlo integral type definition.

	    typedef MC_integral<value_type> integral_type;

	    /* Public constructors: */
	    /*----------------------*/

	    /// Default constructor.

	    MC_integrator():n_calls(0),max_w(0),max_w_eps(0),up_to_date(true),w_hist(NULL),wsum(0),w2sum(0),w3sum(0),w4sum(0),eps(0){}

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    virtual ~MC_integrator()
	    {
		if(w_hist!=NULL)
		{
		    delete w_hist;
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/
	    
	    /// Raises the number of calls by 1 and the weight sums by the
	    /// current weight, without calling the generation method.
	    
	    virtual void refresh_cross_section(bool with_integrand=true)
	    {
		if(!this->validate_weight())
		{
		    return;
		}
		++n_calls;
		if(!with_integrand)
		{
		    update_by_weight(this->weight());
		}
		else
		{
		    if(this->validate_integrand())
		    {
			update_by_weight(weighted_integrand());
		    }
		}
		up_to_date=false;
	    }

	    /// Resets the the state of the MC generator. Resets adaptive
	    /// channels weights and grids as well.

	    virtual void reset()
	    {
		this->weight()=(value_type)0;
		wsum =(value_type)0;
		w2sum=(value_type)0;
		w3sum=(value_type)0;
		w4sum=(value_type)0;
		n_calls=0;
		integral.value=(value_type)0;
		integral.error=std::numeric_limits<value_type>::infinity();
		integral.error_error=std::numeric_limits<value_type>::infinity();
		up_to_date=true;
	    }

	    /// Resets cross section and weight sums of the generator, but keeps
	    /// channel weights and grids intact.

	    virtual void reset_cross_section()
	    {
		this->weight()=(value_type)0;
		wsum =(value_type)0;
		w2sum=(value_type)0;
		w3sum=(value_type)0;
		w4sum=(value_type)0;
		n_calls=0;
		integral.value=(value_type)0;
		integral.error=std::numeric_limits<value_type>::infinity();
		integral.error_error=std::numeric_limits<value_type>::infinity();
		up_to_date=true;
	    }

	    /// Sets up the weight histogram (default 1000 bins).

	    virtual const weight_histogrammer<value_t>* bin_weights(size_type bins)
	    {
		if(w_hist!=NULL)
		{
		    delete w_hist;
		}
		if(max_w==(value_type)0)
		{
		    w_hist=new weight_histogrammer<value_type>(bins);
		}
		else
		{
		    w_hist=new weight_histogrammer<value_type>(max_w,bins);
		}
		return w_hist;
	    }

	    /// Sets the epsilon-reduction to the maximum weight.

	    virtual value_type reduce_max_weight(const value_type& eps_)
	    {
		if(w_hist==NULL)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"cannot reduce max weight--no histogram built"<<endlog;
		    max_w_eps=max_w;
		}
		else
		{
		    eps=eps_;
		    max_w_eps=w_hist->max_weight(eps);
		}
		return max_w_eps;
	    }

	    /// Generates an unweighted event. Internal parameters and the total
	    /// cross section are updated during the rejection loop.

	    template<class rng_t>void generate_unweighted()
	    {
		if(!this->validate_integrand())
		{
		    return;
		}
		value_type rho,x;
		do
		{
		    this->generate();
		    x=weighted_integrand();
		    if(!is_finite_number(x))
		    {
			continue;
		    }
		    this->update();
		    refresh_cross_section();
		    rho=random_number_stream<value_type,rng_t>::throw_number(0,max_w_eps);
		}
		while(rho>x);
		refresh_max_weight();
		this->weight()=this->cross_section().value;
		this->integrand()=(value_type)1;
	    }

	    /// Generates an unweighted event, where the intergrand is
	    /// provided by the argument function object. Internal parameters
	    /// and the cross sectioin are updated during the rejection loop.

	    template<class rng_t,class func_type>void generate_unweighted(const func_type& func)
	    {
		value_type rho,x;
		do
		{
		    this->generate();
		    this->integrand()=func();
		    x=weighted_integrand();
		    if(!is_finite_number(x))
		    {
			continue;
		    }
		    this->update();
		    refresh_cross_section();
		    rho=random_number_stream<value_type,rng_t>::throw_number(0,max_w_eps);
		}
		while(rho>x);
		refresh_max_weight();
		this->weight()=this->cross_section().value;
		this->integrand()=(value_type)1;
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Returns the product of MC weight and integrand.

	    value_type weighted_integrand() const
	    {
		return this->weight()*this->integrand();
	    }

	    /// Returns the maximal received weight x integrand.

	    const value_type& max_weight() const
	    {
		return max_w;
	    }

	    /// Returns the unweighting maximal weight.

	    const value_type& max_weight_epsilon() const
	    {
		return max_w_eps;
	    }

	    /// Cross section computation method. If not up to date, evaluates the mean, variance
	    /// and fourth estimator from the weight sums.

	    virtual MC_integral<value_t> cross_section() const
	    {
		if(up_to_date)
		{
		    return integral;
		}
		value_type N=(value_type)n_calls;
		value_type Ninv=(value_type)1/N;
		value_type result=(n_calls!=0)?Ninv*wsum:(value_type)0;
		MC_integral<value_type>sigma(result);
		value_type N2=N*(N-(value_type)1);
		value_type E2=(w2sum-Ninv*wsum*wsum)/N2;
		sigma.error=(n_calls>1)?std::sqrt(E2):(std::numeric_limits<value_type>::infinity());
		value_type N4=N2*(N-(value_type)2)*(N-(value_type)3);
		value_type N24=N2/N4;
		sigma.error_error=(n_calls>3)?std::pow(N24*Ninv*Ninv*(w4sum-(value_type)4*Ninv*w3sum*wsum+(value_type)3*Ninv*w2sum*w2sum)+((value_type)1-N2*N24)*E2*E2,0.25):(std::numeric_limits<value_type>::infinity());
		integral=sigma;
		up_to_date=true;
		return sigma;
	    }

	    /// Monte Carlo efficiency (in percentage).

	    value_type efficiency() const
	    {
		if(n_calls*max_w==(value_type)0)
		{
		    return 0;
		}
		return (value_type)100*wsum/(n_calls*max_w);
	    }

	    /// Returns the number of calls to the generator.

	    const size_type& calls() const
	    {
		return n_calls;
	    }

	    /// Returns the maximal weight accuracy.

	    const value_type& epsilon() const
	    {
		return eps;
	    }
	    
	    /// Weight histogram plot output. Returns NULL if no weight
	    /// histogram was built.
	    
	    //TODO: get this out
	    plot_script* plot_weights(const std::string& filename,const char* term=NULL) const
	    {
		if(w_hist==NULL)
		{
		    return NULL;
		}
		return w_hist->plot(filename,term);
	    }

	protected:

	    /// Number of cross section update calls.

	    size_type n_calls;

	    /// Monte-Carlo integral.
	    
	    mutable MC_integral<value_type> integral;

	    /// Maximal weight.

	    value_type max_w;

	    /// Epsilon-reduced maximal weight.

	    value_type max_w_eps;

	    /// Time stamp on the Monte-Carlo integral.

	    mutable bool up_to_date;

	    /// Weight histogramming instance.

	    weight_histogrammer<value_t>* w_hist;

	    /// Refreshes the epsilon-reduced maximal weight.

	    void refresh_max_weight()
	    {
		max_w_eps=(w_hist==NULL)?(max_w):(w_hist->max_weight(eps));
	    }

	private:

	    /* Sum of event weights: */

	    value_type wsum;

	    /* Sum of squares of event weights: */
	    
	    value_type w2sum;

	    /* Sum of cubes of event weights: */

	    value_type w3sum;

	    /* Sum of fourth powers of event weights: */

	    value_type w4sum;

	    /* Maximal weight accuracy parameter: */

	    value_type eps;

	    /* Utility helper: */

	    void update_by_weight(const value_type& w)
	    {
		max_w=std::max(w,max_w);
		wsum+=w;
		value_type w2=w*w;
		w2sum+=w2;
		w3sum+=(w2*w);
		w4sum+=(w2*w2);
		if(w_hist!=NULL)
		{
		    w_hist->insert(w);
		}
	    }
    };

    /// Wrapper class, converting a generator into an integrator.

    template<class value_t>class MC_generator_wrapper: public MC_integrator<value_t>
    {
	public:

	    typedef value_t value_type;
	    typedef std::size_t size_type;
	    typedef MC_generator<value_t> generator_type;
	    typedef MC_integrator<value_t> base_type;

	    /// Constructor.

	    MC_generator_wrapper(generator_type* generator_,bool owner_=true):generator(generator_),integrator(dynamic_cast<MC_integrator<value_t>*>(generator_)),owner(owner_){}

	    /// Destructor.

	    ~MC_generator_wrapper()
	    {
		if(owner)
		{
		    delete generator;
		}
	    }

	    bool generate()
	    {
		bool result=generator->generate();
		this->weight()=generator->weight();
		return result;
	    }

	    bool evaluate_weight()
	    {
		bool result=generator->evaluate_weight();
		this->weight()=generator->weight();
		return result;
	    }

	    void update()
	    {
		if(integrator!=NULL)
		{
		    integrator->integrand()=this->integrand();
		    integrator->update();
		}
	    }

	    void adapt()
	    {
		if(integrator!=NULL)
		{
		    integrator->adapt();
		}
	    }

	private:

	    MC_generator<value_t>* const generator;
	    MC_integrator_base<value_t>* const integrator;
	    bool owner;
    };
}

#endif /*CAMGEN_MC_GEN_H_*/

