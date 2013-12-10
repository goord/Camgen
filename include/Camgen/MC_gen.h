//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file MC_gen.h
    \brief Abstract base class template for Monte Carlo generator classes.
 */

#ifndef CAMGEN_MC_GEN_H_
#define CAMGEN_MC_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for the Monte Carlo generators. *
 *                                                     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <limits>
#include <Camgen/utils.h>
#include <Camgen/logstream.h>
#include <Camgen/MC_int.h>
#include <Camgen/whist.h>
#include <Camgen/rn_strm.h>

namespace Camgen
{
    template<class value_t>class MC_generator
    {
	public:

	    /* Public type definitions: */
	    /*--------------------------*/

	    typedef std::size_t size_type;

	    /// Value type definition.
	    
	    typedef value_t value_type;

	    /// Evaluated Monte Carlo integral type definition.

	    typedef MC_integral<value_type> integral_type;

	    /* Public static methods: */
	    /*------------------------*/

	    /// Static function checking whether the argument is defined and not
	    /// infinite.

	    static bool valid(const value_type& x)
	    {
		return (x==x and x!=std::numeric_limits<value_type>::infinity() and x!=-std::numeric_limits<value_type>::infinity());
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /// Default constructor.

	    MC_generator():n_calls(0),max_w(0),max_w_eps(0),up_to_date(true),w_hist(NULL),w(0),wsum(0),w2sum(0),w3sum(0),w4sum(0),f(1),eps(0){}

	    /* Destructor: */
	    /*-------------*/

	    /// Destructor.

	    virtual ~MC_generator()
	    {
		if(w_hist!=NULL)
		{
		    delete w_hist;
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/
	    
	    /// 'Bare generation' method (purely virtual). Should update the weight
	    /// data member.
	    
	    virtual bool generate()=0;

	    /// Event weight calculation method. Trivial by default.

	    virtual bool evaluate_weight()
	    {
		return true;
	    }

	    /// Update method. Should refresh internal parameters with
	    /// information given by integrand.

	    virtual void update(){}

	    /// Adaptation method. Refines internal grids to improve efficiency. Trivial by default.

	    virtual void adapt(){}

	    /// Operator combining the generation cross section refreshing and
	    /// updating methods.

	    value_type operator()()
	    {
		generate();
		refresh_cross_section();
		update();
		return w*f;
	    }

	    /// Raises the number of calls by 1 and the weight sums by the
	    /// current weight, without calling the generation method.
	    
	    virtual void refresh_cross_section(bool with_integrand=true)
	    {
		if(!valid(w))
		{
		    return;
		}
		++n_calls;
		if(!with_integrand)
		{
		    update_by_weight(w);
		}
		else
		{
		    if(valid(this->integrand()))
		    {
			update_by_weight(w*f);
		    }
		}
		up_to_date=false;
	    }

	    /// Resets the the state of the MC generator. Resets adaptive
	    /// channels weights and grids as well.

	    virtual void reset()
	    {
		w=(value_type)0;
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
		w=(value_type)0;
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
		if(!valid(f))
		{
		    return;
		}
		value_type rho,x;
		do
		{
		    generate();
		    x=w*f;
		    if(!valid(x))
		    {
			continue;
		    }
		    update();
		    refresh_cross_section();
		    rho=random_number_stream<value_type,rng_t>::throw_number(0,max_w_eps);
		}
		while(rho>x);
		refresh_max_weight();
		w=this->cross_section().value;
		f=(value_type)1;
	    }

	    /// Generates an unweighted event, where the intergrand is
	    /// provided by the argument function object. Internal parameters
	    /// and the cross sectioin are updated during the rejection loop.

	    template<class rng_t,class func_type>void generate_unweighted(const func_type& func)
	    {
		value_type rho,x;
		do
		{
		    generate();
		    f=func();
		    x=w*f;
		    if(!valid(x))
		    {
			continue;
		    }
		    update();
		    refresh_cross_section();
		    rho=random_number_stream<value_type,rng_t>::throw_number(0,max_w_eps);
		}
		while(rho>x);
		refresh_max_weight();
		w=this->cross_section().value;
		f=(value_type)1;
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Returns the Monte Carlo weight reference.

	    value_type& weight()
	    {
		return w;
	    }

	    /// Returns the Monte Carlo const weight reference.

	    const value_type& weight() const
	    {
		return w;
	    }

	    /// Returns the Monte Carlo weight address.

	    value_type* get_weight()
	    {
		return &w;
	    }

	    /// Returns the Monte Carlo const weight reference.

	    const value_type* get_weight() const
	    {
		return &w;
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

	    /// Returns the reference to the integrand (1 by default).

	    value_type& integrand()
	    {
		return f;
	    }

	    /// Returns the const reference to the integrand (1 by default).

	    const value_type& integrand() const
	    {
		return f;
	    }

	    /// Returns the integrand address.

	    value_type* get_integrand()
	    {
		return &f;
	    }

	    /// Returns the const pointer to the integrand.

	    const value_type* get_integrand() const
	    {
		return &f;
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

	    /* Serialisation: */
	    /*----------------*/

	    /// Loads cross section data from input stream.

	    virtual std::istream& load(std::istream& is)
	    {
		std::string initflag;
		do
		{
		    std::getline(is,initflag);
		}
		while(initflag!="<init>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read"<<endlog;
		    return is;
		}
		safe_read(is,integral.value);
		safe_read(is,integral.error);
		safe_read(is,integral.error_error);
		is>>n_calls;
		safe_read(is,wsum);
		safe_read(is,w2sum);
		safe_read(is,w3sum);
		safe_read(is,w4sum);
		safe_read(is,max_w);
		safe_read(is,max_w_eps);
		safe_read(is,eps);
		is>>std::ws;
		if(is.peek()=='N')
		{
		    if(w_hist!=NULL)
		    {
			delete w_hist;
			w_hist=NULL;
		    }
		    is.ignore(1);
		}
		else
		{
		    if(w_hist!=NULL)
		    {
			delete w_hist;
		    }
		    w_hist=weight_histogrammer<value_t>::create_instance(is);
		}
		std::string endflag;
		do
		{
		    std::getline(is,endflag);
		}
		while(endflag!="</init>" and !is.eof());
		return is;
	    }

	    /// Saves cross section data from input stream.

	    virtual std::ostream& save(std::ostream& os) const
	    {
		os<<"<init>"<<std::endl;
		safe_write(os,integral.value);
		os<<"\t";
		safe_write(os,integral.error);
		os<<"\t";
		safe_write(os,integral.error_error);
		os<<std::endl;
		os<<n_calls;
		os<<"\t";
		safe_write(os,wsum);
		os<<"\t";
		safe_write(os,w2sum);
		os<<"\t";
		safe_write(os,w3sum);
		os<<"\t";
		safe_write(os,w4sum);
		os<<"\t";
		safe_write(os,max_w);
		os<<"\t";
		safe_write(os,max_w_eps);
		os<<"\t";
		safe_write(os,eps);
		os<<std::endl;
		if(w_hist==NULL)
		{
		    os<<'N'<<std::endl;
		}
		else
		{
		    w_hist->save(os);
		}
		os<<"</init>"<<std::endl;
		return os;
	    }
	    
	    /// Weight histogram plot output. Returns NULL if no weight
	    /// histogram was built.
	    
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

	    /* Event weight: */

	    value_type w;

	    /* Sum of event weights: */

	    value_type wsum;

	    /* Sum of squares of event weights: */
	    
	    value_type w2sum;

	    /* Sum of cubes of event weights: */

	    value_type w3sum;

	    /* Sum of fourth powers of event weights: */

	    value_type w4sum;

	    /* Integrand: */

	    value_type f;

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
}

#endif /*CAMGEN_MC_GEN_H_*/

