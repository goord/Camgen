//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_HAD_IS_H_
#define CAMGEN_HAD_IS_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Hadronic initial state generator class using the LHAPDF weights,        *
 * combined with the parni adaptive grid to optimize x0,x1 generation with *
 * respect to some integrand.                                              *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/Minkowski.h>
#include <Camgen/init_state.h>
#include <Camgen/pdf_wrapper.h>
#include <Camgen/rn_strm.h>
#include <Camgen/power_law.h>
#include <Camgen/inv_cosh.h>
#include <Camgen/ps_vol.h>

namespace Camgen
{
    /* LHAPDF-weighted adaptive (anti-)proton-(anti-)proton initial state
    generator, parameterised by the momentum fractions carried by the
    initial partons. If the LHAPDF headers were not included, the hadron
    structure functions are set to unity. */

    template<class model_t,class rng_t>class hadronic_is_xx<model_t,rng_t,Minkowski_type>: public initial_state<model_t,2>
    {
	typedef initial_state<model_t,2> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::size_type size_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Static integer constants: */
	    
	    static const size_type k0=spacetime_type::timelike_direction;
	    static const size_type kL=(model_type::beam_direction<0)?(-model_type::beam_direction):(model_type::beam_direction);

	    /* Public data members */
	    /*---------------------*/

	    /* First parton flavour: */

	    const int hadron1;

	    /* Second parton flavour: */

	    const int hadron2;

	    /* Public constructors */
	    /*---------------------*/

	    /* Constructor with anti-proton flags (true denotes pbar-beam): */

	    hadronic_is_xx(bool q1,bool q2):base_type(true,true,true),hadron1(q1?(-2212):2212),hadron2(q2?(-2212):2212),flip(false)
	    {
		pdf_wrapper::initialise(MC_config::pdf_name(),MC_config::pdf_number());
		xmin.assign(static_cast<value_type>(pdf_wrapper::xmin()));
		xmax.assign(static_cast<value_type>(pdf_wrapper::xmax()));
		xgen=new parni<value_type,2,rng_t>(&xvec,xmin,xmax,2*grid_bins(),grid_mode());
	    }

	    /* Copy constructor: */
	    
	    hadronic_is_xx(const hadronic_is_xx<model_t,rng_t>& other):base_type(other),hadron1(other.hadron1),hadron2(other.hadron2),flip(false)
	    {
		pdf_wrapper::initialise(MC_config::pdf_name(),MC_config::pdf_number());
		xmin.assign(static_cast<value_type>(pdf_wrapper::xmin()));
		xmax.assign(static_cast<value_type>(pdf_wrapper::xmax()));
		xgen=new parni<value_type,2,rng_t>(&xvec,xmin,xmax,2*grid_bins(),grid_mode());
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /* Destructor: */

	    ~hadronic_is_xx()
	    {
		delete xgen;
	    }

	    /* Public modifiers */
	    /*------------------*/

	    /* Hadronic invariant mass refresher (incoming hadrons are assumed
	     * to be massless). */

	    bool refresh_Ecm()
	    {
		value_type E1=this->beam_energy(0);
		value_type E2=this->beam_energy(1);
		return this->set_s((value_type)4*E1*E2);
	    }

	    /* Generation method implementation. */
	    
	    bool generate()
	    {
		if(!xgen->generate())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		
		value_type P1=xvec[0]*this->beam_energy(0);
		value_type P2=xvec[1]*this->beam_energy(1);
		
		flip=rn_stream::throw_coin();
		size_type pa=flip?0:1;
		size_type pb=flip?1:0;
		
		value_type Ehat1=std::sqrt(P1*P1+this->m2(pa));
		value_type Ehat2=std::sqrt(P2*P2+this->m2(pb));
		if(!this->set_s_hat(this->m2(0)+this->m2(1)+(value_type)2*(Ehat1*Ehat2+P1*P2)))
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}

		this->p(pa).assign(0);
		this->p(pa)[k0]=Ehat1;
		this->p(pa)[kL]=P1;

		this->p(pb).assign(0);
		this->p(pb)[k0]=Ehat2;
		this->p(pb)[kL]=-P2;

		this->weight()=xgen->weight();
		return true;
	    }

	    /* Weight evaluation implementation: */
	    
	    bool evaluate_weight()
	    {
		if(!xgen->evaluate_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		const momentum_type& p1=this->p(0);
		const momentum_type& p2=this->p(1);
		if(p1[0]<this->m(0) or p2[0]<this->m(1))
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		value_type P1=std::abs(p1[kL]);
		value_type P2=std::abs(p2[kL]);
		if(!this->set_s_hat(this->m2(0)+this->m2(1)+(value_type)2*(p1[0]*p2[0]+P1*P2)))
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		xvec[0]=P1/this->beam_energy(0);
		xvec[1]=P2/this->beam_energy(1);
		if(!xgen->evaluate_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->weight()=xgen->weight();
		return true;
	    }

	    /* Flux factor output: */

	    value_type flux_factor() const
	    {
		size_type pa=flip?0:1;
		size_type pb=flip?1:0;
		int q1=(this->parton_id(pa)==21)?0:((hadron1<0)?(-this->parton_id(pa)):(this->parton_id(pa)));
		int q2=(this->parton_id(pb)==21)?0:((hadron2<0)?(-this->parton_id(pb)):(this->parton_id(pb)));

		return static_cast<value_type>(pdf_wrapper::ff(xvec[0],xvec[1],q1,q2,this->mu_F))/std::sqrt(Kallen(this->s_hat(),this->m2(0),this->m2(1)));
	    }

	    /* Update method: */

	    void update()
	    {
		xgen->integrand()=this->integrand()/xgen->weight();
		xgen->update();
	    }
	    
	    /* Adaptor method implementation: */

	    void adapt_grids()
	    {
		xgen->adapt();
	    }

	    /* Resets adaptive grids: */

	    void reset()
	    {
		this->base_type::reset();
		xgen->reset();
	    }

	    /* Resets adaptive grids: */

	    void reset_cross_section()
	    {
		this->base_type::reset_cross_section();
		xgen->reset_cross_section();
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /* Clone method implementation: */

	    hadronic_is_xx<model_t,rng_t,Minkowski_type>* clone() const
	    {
		return new hadronic_is_xx<model_t,rng_t,Minkowski_type>(*this);
	    }

	    /* Beam particle id: */

	    int beam_id(size_type i) const
	    {
		if(i<2)
		{
		    return (i==0)?hadron1:hadron2;
		}
		return 0;
	    }

	    /* Cernlib pdf group number: */

	    int pdfg(size_type i) const
	    {
		return pdf_wrapper::group();
	    }

	    /* Cernlib pdf set number: */

	    int pdfs(size_type i) const
	    {
		return pdf_wrapper::set();
	    }

	    /* Returns the coupling constant corresponding to the factorization
	    scale: */

	    value_type alpha_s() const
	    {
		return static_cast<value_type>(pdf_wrapper::alpha_s(this->mu_F));
	    }

	    /* x-value redout: */

	    const value_type& x(size_type i) const
	    {
		return xvec[i];
	    }

	    /* Lower bounds readout: */

	    const value_type& x_min(size_type i) const
	    {
		return xmin[i];
	    }

	    /* Upper bound readout: */

	    const value_type& x_max(size_type i) const
	    {
		return xmax[i];
	    }

	    /* Serialisation: */
	    /*----------------*/

	    /* Returns the polymorphic type. */

	    std::string type() const
	    {
		switch(hadron1)
		{
		    case 2212:
			if(hadron2==2212)
			{
			    return "pp_xx";
			}
			else if(hadron2==-2212)
			{
			    return "ppbar_xx";
			}
			else
			{
			    return "X";
			}
			break;
		    case -2212:
			if(hadron2==2212)
			{
			    return "pbarp_xx";
			}
			else if(hadron2==-2212)
			{
			    return "pbarpbar_xx";
			}
			else
			{
			    return "X";
			}
			break;
		    default:
			return "X";
		}
	    }

	    /// Overloaded settings printing function.

	    std::ostream& print_settings(std::ostream& os) const
	    {
		os<<std::setw(30)<<std::left<<"beam energy -1:"<<this->beam_energy(0)<<std::endl;
		os<<std::setw(30)<<std::left<<"hadron -1:"<<hadron1<<std::endl;
		os<<std::setw(30)<<std::left<<"beam energy -2:"<<this->beam_energy(0)<<std::endl;
		os<<std::setw(30)<<std::left<<"hadron -2:"<<hadron2<<std::endl;
		os<<std::setw(30)<<std::left<<"Ecm:"<<this->Ecm()<<std::endl;
		os<<std::setw(30)<<std::left<<"PDF set:"<<pdf_name()<<std::endl;
		os<<std::setw(30)<<std::left<<"PDF nr:"<<pdf_number()<<std::endl;
		os<<std::setw(30)<<std::left<<"mu_F:"<<this->mu_F<<std::endl;
		os<<std::setw(30)<<std::left<<"alpha_s:"<<alpha_s()<<std::endl;
		os<<std::setw(30)<<std::left<<"mhatmin:"<<this->m_hat_min()<<std::endl;
		os<<std::setw(30)<<std::left<<"x1 interval:"<<"["<<xmin[0]<<','<<xmax[0]<<']'<<std::endl;
		os<<std::setw(30)<<std::left<<"x2 interval:"<<"["<<xmin[1]<<','<<xmax[1]<<']'<<std::endl;
		return os;
	    }

	    /* Derived loading method: */

	    std::istream& load_data(std::istream& is)
	    {
		safe_read(is,xmin[0]);
		safe_read(is,xmin[1]);
		safe_read(is,xmax[0]);
		safe_read(is,xmax[1]);
		if(xgen!=NULL)
		{
		    delete xgen;
		}
		xgen=parni<value_type,2,rng_t>::create_instance(&xvec,is);
		return is;
	    }

	    /* Derived saving method: */

	    std::ostream& save_data(std::ostream& os) const
	    {
		safe_write(os,xmin[0]);
		os<<"\t";
		safe_write(os,xmin[1]);
		os<<"\t\t";
		safe_write(os,xmax[0]);
		os<<"\t";
		safe_write(os,xmax[1]);
		os<<std::endl;
		xgen->save(os);
		return os;
	    }

	private:

	    /* x-values, minima and maxima: */
	    
	    vector<value_type,2>xvec,xmin,xmax;
	    
	    /* Adaptive grid generating x-values: */
	    
	    parni<value_type,2,rng_t>* xgen;

	    /* Boolean denoting whether the initial state was swapped: */

	    bool flip;
    };
    template<class model_t,class rng_t>const typename hadronic_is_xx<model_t,rng_t,Minkowski_type>::size_type hadronic_is_xx<model_t,rng_t,Minkowski_type>::k0;
    template<class model_t,class rng_t>const typename hadronic_is_xx<model_t,rng_t,Minkowski_type>::size_type hadronic_is_xx<model_t,rng_t,Minkowski_type>::kL;
    

    /* LHAPDF-weighted adaptive (anti-)proton-(anti-)proton initial state
    generator, parameterized by the partonic invariant mass and the rapidity
    of the partonic system. If the LHAPDF-header was not included at
    configuration, the hadronic structure function are set to unity. */

    template<class model_t,class rng_t>class hadronic_is_sy<model_t,rng_t,Minkowski_type>: public initial_state<model_t,2>
    {
	/* Type definitions: */

	typedef initial_state<model_t,2> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::size_type size_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Static integer constants: */
	    
	    static const size_type k0=spacetime_type::timelike_direction;
	    static const size_type kL=(model_type::beam_direction<0)?(-model_type::beam_direction):(model_type::beam_direction);

	    /* Public data members */
	    /*---------------------*/

	    /* First parton flavour: */

	    const int hadron1;

	    /* Second beam hadron flavour: */

	    const int hadron2;

	    /// Exponent in power-law invariant mass generation.
	    
	    value_type nu_tau;

	    /* Public constructors */
	    /*---------------------*/

	    /* Constructor with anti-proton flags (true denotes pbar-beam): */

	    hadronic_is_sy(bool q1,bool q2):base_type(true,true,true),hadron1(q1?(-2212):2212),hadron2(q2?(-2212):2212),nu_tau(shat_exponent()),flip(false)
	    {
		pdf_wrapper::initialise(MC_config::pdf_name(),MC_config::pdf_number());
		LHA_ymax=static_cast<value_type>(pdf_wrapper::ymax());
		ymax=LHA_ymax;
		tau_gen=new adaptive_s_generator<value_type,rng_t>(&tau,new pl_s_generator<value_type,rng_t>(NULL,&nu_tau));
		tau_gen->set_s_min_min(0);
		tau_gen->set_s_min(0);
		tau_gen->set_s_max_max(1);
		tau_gen->set_s_max(1);
		y_gen=new adaptive_s_generator<value_type,rng_t>(&y,new inv_cosh_y_generator<value_type,rng_t>);
		if(ymax!=std::numeric_limits<value_type>::infinity())
		{
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		}
	    }

	    /* Copy constructor: */
	    
	    hadronic_is_sy(const hadronic_is_sy<model_t,rng_t>& other):base_type(other),hadron1(other.hadron1),hadron2(other.hadron2),nu_tau(other.nu_tau),flip(false)
	    {
		pdf_wrapper::initialise(MC_config::pdf_name(),MC_config::pdf_number());
		LHA_ymax=static_cast<value_type>(pdf_wrapper::ymax());
		ymax=LHA_ymax;
		tau_gen=new adaptive_s_generator<value_type,rng_t>(&tau,new pl_s_generator<value_type,rng_t>(NULL,&nu_tau));
		tau_gen->set_s_min_min(0);
		tau_gen->set_s_min(0);
		tau_gen->set_s_max_max(1);
		tau_gen->set_s_max(1);
		y_gen=new adaptive_s_generator<value_type,rng_t>(&y,new inv_cosh_y_generator<value_type,rng_t>);
		if(LHA_ymax!=std::numeric_limits<value_type>::infinity())
		{
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		}
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /* Destructor: */

	    ~hadronic_is_sy()
	    {
		delete tau_gen;
		delete y_gen;
	    }

	    /* Public modifiers */
	    /*------------------*/
	    
	    /* Hadronic invariant mass refresher (incoming hadrons are assumed
	     * to be massless). */

	    bool refresh_Ecm()
	    {
		value_type E1=this->beam_energy(0);
		value_type E2=this->beam_energy(1);
		if(this->set_s((value_type)4*E1*E2))
		{
		    tau_gen->set_s_min_min((this->s_hat_min()-this->m2(0)-this->m2(1))/this->s());
		    tau_gen->set_s_min((this->s_hat_min()-this->m2(0)-this->m2(1))/this->s());
		    ymax=std::min(-(value_type)0.5*std::log(tau_gen->s_min()),LHA_ymax);
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		    return true;
		}
		return false;
	    }

	    /* Sets the minimal partonic invariant mass-squared to the argument. */
	    
	    bool set_s_hat_min(const value_type& s)
	    {
		if(this->base_type::set_s_hat_min(s))
		{
		    tau_gen->set_s_min_min((this->s_hat_min()-this->m2(0)-this->m2(1))/this->s());
		    tau_gen->set_s_min((this->s_hat_min()-this->m2(0)-this->m2(1))/this->s());
		    ymax=std::min(-(value_type)0.5*std::log(tau_gen->s_min()),LHA_ymax);
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		    return true;
		}
		return false;
	    }

	    /* Sets the minimal partonic invariant mass to the argument. */
	    
	    bool set_m_hat_min(const value_type& m)
	    {
		if(this->base_type::set_m_hat_min(m))
		{
		    tau_gen->set_s_min_min((this->s_hat_min()-this->m2(0)-this->m2(1))/this->s());
		    tau_gen->set_s_min((this->s_hat_min()-this->m2(0)-this->m2(1))/this->s());
		    ymax=std::min(-(value_type)0.5*std::log(tau_gen->s_min()),LHA_ymax);
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		    return true;
		}
		return false;
	    }

	    /* Generation method implementation. */
	    
	    bool generate()
	    {
		if(!tau_gen->generate())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"tau generation failed within ["<<tau_gen->s_min()<<','<<tau_gen->s_max()<<']'<<endlog;
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type ybound=std::min(ymax,-(value_type)0.5*std::log(tau));
		y_gen->set_s_range(-ybound,ybound);
		if(!y_gen->generate())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"y generation failed within ["<<y_gen->s_min()<<','<<y_gen->s_max()<<']'<<endlog;
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type sqrttau=std::sqrt(tau);
		value_type expy=std::exp(y);
		x1=sqrttau*expy;
		x2=sqrttau/expy;
		
		value_type P1=x1*this->beam_energy(0);
		value_type P2=x2*this->beam_energy(1);
		
		flip=rn_stream::throw_coin();
		size_type pa=flip?0:1;
		size_type pb=flip?1:0;
		
		value_type Ehat1=std::sqrt(P1*P1+this->m2(pa));
		value_type Ehat2=std::sqrt(P2*P2+this->m2(pb));

		this->p(pa).assign(0);
		this->p(pa)[k0]=Ehat1;
		this->p(pa)[kL]=P1;
		this->p(pb).assign(0);
		this->p(pb)[k0]=Ehat2;
		this->p(pb)[kL]=-P2;

		if(!this->set_s_hat(this->m2(0)+this->m2(1)+(value_type)2*(Ehat1*Ehat2+P1*P2)))
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		this->weight()=tau_gen->weight()*y_gen->weight();
		return true;
	    }
	    
	    /* Weight evaluation implementation: */
	    
	    bool evaluate_weight()
	    {
		const momentum_type& p1=this->p(0);
		const momentum_type& p2=this->p(1);
		if(p1[k0]<this->m(0) or p2[k0]<this->m(1))
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		value_type P1=std::abs(p1[kL]);
		value_type P2=std::abs(p2[kL]);
		if(!this->set_s_hat(this->m2(0)+this->m2(1)+(value_type)2*(p1[k0]*p2[k0]+P1*P2)))
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		x1=P1/this->beam_energy(0);
		x2=P2/this->beam_energy(1);
		tau=x1*x2;
		if(!tau_gen->evaluate_weight())
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		y=(value_type)0.5*std::log(x1/x2);
		value_type ybound=std::max(ymax,-(value_type)0.5*std::log(tau));
		y_gen->set_s_range(-ybound,ybound);
		if(!y_gen->evaluate_weight())
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		this->weight()=tau_gen->weight()*y_gen->weight();
		return true;
	    }

	    /* Flux factor output: */

	    value_type flux_factor() const
	    {
		size_type pa=flip?0:1;
		size_type pb=flip?1:0;
		int q1=(this->parton_id(pa)==21)?0:((hadron1<0)?(-this->parton_id(pa)):(this->parton_id(pa)));
		int q2=(this->parton_id(pb)==21)?0:((hadron2<0)?(-this->parton_id(pb)):(this->parton_id(pb)));
		
		return static_cast<value_type>(pdf_wrapper::ff(x1,x2,q1,q2,this->mu_F))*x1*x2/(tau*std::sqrt(Kallen(this->s_hat(),this->m2(0),this->m2(1))));
	    }

	    /* Update method: */

	    void update()
	    {
		tau_gen->integrand()=this->integrand();
		tau_gen->update_weight();
		y_gen->integrand()=this->integrand();
		y_gen->update_weight();
	    }
	    
	    /* Adaptor method implementation: */

	    void adapt_grids()
	    {
		tau_gen->adapt();
		y_gen->adapt();
	    }

	    /* Resets adaptive grids: */

	    void reset()
	    {
		this->base_type::reset();
		tau_gen->reset();
		y_gen->reset();
	    }

	    /* Resets cross sections: */

	    void reset_cross_section()
	    {
		this->base_type::reset_cross_section();
		tau_gen->reset_cross_section();
		y_gen->reset_cross_section();
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /* Clone method implementation: */

	    hadronic_is_sy<model_t,rng_t,Minkowski_type>* clone() const
	    {
		return new hadronic_is_sy<model_t,rng_t,Minkowski_type>(*this);
	    }

	    /* Beam particle id: */

	    int beam_id(size_type i) const
	    {
		if(i<2)
		{
		    return (i==0)?hadron1:hadron2;
		}
		return 0;
	    }

	    /* Cernlib pdf group number: */

	    int pdfg(size_type i) const
	    {
		return pdf_wrapper::group();
	    }

	    /* Cernlib pdf set number: */

	    int pdfs(size_type i) const
	    {
		return pdf_wrapper::set();
	    }

	    /* Returns the coupling constant corresponding to the factorization
	    scale: */

	    value_type alpha_s() const
	    {
		return pdf_wrapper::alpha_s(this->mu_F);
	    }

	    /* Serialisation: */
	    /*----------------*/

	    /* Returns the polymorphic type. */

	    std::string type() const
	    {
		switch(hadron1)
		{
		    case 2212:
			if(hadron2==2212)
			{
			    return "pp_sy";
			}
			else if(hadron2==-2212)
			{
			    return "ppbar_sy";
			}
			else
			{
			    return "X";
			}
			break;
		    case -2212:
			if(hadron2==2212)
			{
			    return "pbarp_sy";
			}
			else if(hadron2==-2212)
			{
			    return "pbarpbar_sy";
			}
			else
			{
			    return "X";
			}
			break;
		    default:
			return "X";
		}
	    }

	    /// Overloaded settings printing function.

	    std::ostream& print_settings(std::ostream& os) const
	    {
		os<<std::setw(30)<<std::left<<"beam energy -1:"<<this->beam_energy(0)<<std::endl;
		os<<std::setw(30)<<std::left<<"hadron -1:"<<hadron1<<std::endl;
		os<<std::setw(30)<<std::left<<"beam energy -2:"<<this->beam_energy(0)<<std::endl;
		os<<std::setw(30)<<std::left<<"hadron -2:"<<hadron2<<std::endl;
		os<<std::setw(30)<<std::left<<"Ecm:"<<this->Ecm()<<std::endl;
		os<<std::setw(30)<<std::left<<"PDF set:"<<pdf_name()<<std::endl;
		os<<std::setw(30)<<std::left<<"PDF nr:"<<pdf_number()<<std::endl;
		os<<std::setw(30)<<std::left<<"mu_F:"<<this->mu_F<<std::endl;
		os<<std::setw(30)<<std::left<<"alpha_s:"<<alpha_s()<<std::endl;
		os<<std::setw(30)<<std::left<<"mhatmin:"<<this->m_hat_min()<<std::endl;
		os<<std::setw(30)<<std::left<<"shat-sampling exp:"<<nu_tau<<std::endl;
		return os;
	    }

	    /* Derived loading method: */

	    std::istream& load_data(std::istream& is)
	    {
		safe_read(is,nu_tau);
		safe_read(is,ymax);
		safe_read(is,LHA_ymax);
		if(tau_gen!=NULL)
		{
		    delete tau_gen;
		}
		tau_gen=s_generator<value_type,rng_t>::create_instance(&tau,is,NULL,NULL,&nu_tau);
		if(y_gen!=NULL)
		{
		    delete y_gen;
		}
		y_gen=s_generator<value_type,rng_t>::create_instance(&y,is,NULL,NULL,NULL);
		return is;
	    }

	    /* Derived saving method: */

	    std::ostream& save_data(std::ostream& os) const
	    {
		safe_write(os,nu_tau);
		os<<"\t";
		safe_write(os,ymax);
		os<<"\t";
		safe_write(os,LHA_ymax);
		os<<std::endl;
		tau_gen->save(os);
		y_gen->save(os);
		return os;
	    }

	private:
	    
	    value_type tau,y,x1,x2,ymax,LHA_ymax;
	    s_generator<value_type,rng_t>* tau_gen;
	    s_generator<value_type,rng_t>* y_gen;
	    bool flip;
    };
    template<class model_t,class rng_t>const typename hadronic_is_sy<model_t,rng_t,Minkowski_type>::size_type hadronic_is_sy<model_t,rng_t,Minkowski_type>::k0;
    template<class model_t,class rng_t>const typename hadronic_is_sy<model_t,rng_t,Minkowski_type>::size_type hadronic_is_sy<model_t,rng_t,Minkowski_type>::kL;
    
    /* LHAPDF-weighted adaptive (anti-)proton-(anti-)proton initial state
    generator, parameterized by the rapidity of the partonic system. The
    partonic invariant mass should be generated by the final state. */
    
    template<class model_t,class rng_t>class hadronic_is_y<model_t,rng_t,Minkowski_type>: public initial_state<model_t,2>
    {
	/* Type definitions: */

	typedef initial_state<model_t,2> base_type;

	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::spacetime_type spacetime_type;
	    typedef typename base_type::size_type size_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_type,rng_t> rn_stream;

	    /* Static integer constants: */
	    
	    static const size_type k0=spacetime_type::timelike_direction;
	    static const size_type kL=(model_type::beam_direction<0)?(-model_type::beam_direction):(model_type::beam_direction);

	    /* Public data members */
	    /*---------------------*/

	    /* First parton flavour: */

	    const int hadron1;

	    /* Second beam hadron flavour: */

	    const int hadron2;

	    /* Public constructors */
	    /*---------------------*/

	    /* Constructor with anti-proton flags (true denotes pbar-beam): */

	    hadronic_is_y(bool q1,bool q2):base_type(true,true,false),hadron1(q1?(-2212):2212),hadron2(q2?(-2212):2212),flip(false)
	    {
		pdf_wrapper::initialise(MC_config::pdf_name(),MC_config::pdf_number());
		value_type xmin=static_cast<value_type>(pdf_wrapper::xmin());
		y_gen=new adaptive_s_generator<value_type,rng_t>(&y,new inv_cosh_y_generator<value_type,rng_t>);
		if(xmin>(value_type)0)
		{
		    LHA_ymax=(value_type)0.5*std::log(static_cast<value_type>(pdf_wrapper::xmax())/xmin);
		    y_gen->set_s_min_min(-LHA_ymax);
		    y_gen->set_s_min(-LHA_ymax);
		    y_gen->set_s_max_max(LHA_ymax);
		    y_gen->set_s_max(LHA_ymax);
		}
		else
		{
		    LHA_ymax=std::numeric_limits<value_type>::infinity();
		}
		ymax=LHA_ymax;
	    }

	    /* Copy constructor: */
	    
	    hadronic_is_y(const hadronic_is_y<model_t,rng_t>& other):base_type(other),hadron1(other.hadron1),hadron2(other.hadron2),flip(false)
	    {
		pdf_wrapper::initialise(MC_config::pdf_name(),MC_config::pdf_number());
		value_type xmin=static_cast<value_type>(pdf_wrapper::xmin());
		y_gen=new adaptive_s_generator<value_type,rng_t>(&y,new inv_cosh_y_generator<value_type,rng_t>);
		if(xmin>(value_type)0)
		{
		    LHA_ymax=(value_type)0.5*std::log(static_cast<value_type>(pdf_wrapper::xmax())/xmin);
		    y_gen->set_s_min_min(-LHA_ymax);
		    y_gen->set_s_min(-LHA_ymax);
		    y_gen->set_s_max_max(LHA_ymax);
		    y_gen->set_s_max(LHA_ymax);
		}
		else
		{
		    LHA_ymax=std::numeric_limits<value_type>::infinity();
		}
		ymax=LHA_ymax;
	    }

	    /* Public destructors */
	    /*--------------------*/

	    /* Destructor: */

	    ~hadronic_is_y()
	    {
		delete y_gen;
	    }

	    /* Public modifiers */
	    /*------------------*/
	    
	    /* Hadronic invariant mass refresher (incoming hadrons are assumed
	     * to be massless). */

	    bool refresh_Ecm()
	    {
		value_type E1=this->beam_energy(0);
		value_type E2=this->beam_energy(1);
		if(this->set_s((value_type)4*E1*E2))
		{
		    tau_min=(this->s_hat_min()-this->m2(0)-this->m2(1))/this->s();
		    ymax=std::min(-(value_type)0.5*std::log(tau_min),LHA_ymax);
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		    return true;
		}
		return false;
	    }

	    /* Sets the minimal partonic invariant mass-squared to the argument. */
	    
	    bool set_s_hat_min(const value_type& s)
	    {
		if(this->base_type::set_s_hat_min(s))
		{
		    tau_min=(this->s_hat_min()-this->m2(0)-this->m2(1))/this->s();
		    ymax=std::min(-(value_type)0.5*std::log(tau_min),LHA_ymax);
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		    return true;
		}
		return false;
	    }

	    /* Sets the minimal partonic invariant mass to the argument. */
	    
	    bool set_m_hat_min(const value_type& m)
	    {
		if(this->base_type::set_m_hat_min(m))
		{
		    tau_min=(this->s_hat_min()-this->m2(0)-this->m2(1))/this->s();
		    ymax=std::min(-(value_type)0.5*std::log(tau_min),LHA_ymax);
		    y_gen->set_s_min_min(-ymax);
		    y_gen->set_s_min(-ymax);
		    y_gen->set_s_max_max(ymax);
		    y_gen->set_s_max(ymax);
		    return true;
		}
		return false;
	    }

	    /* Generation method implementation. */
	    
	    bool generate()
	    {
		tau=this->s_hat()/this->s();
		value_type ybound=std::min(ymax,-(value_type)0.5*std::log(tau));
		y_gen->set_s_min(-ybound);
		y_gen->set_s_max(ybound);
		if(!y_gen->generate())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		value_type sqrttau=std::sqrt(tau);
		value_type expy=std::exp(y);
		x1=sqrttau*expy;
		x2=sqrttau/expy;
		value_type P1=x1*this->beam_energy(0);
		value_type P2=x2*this->beam_energy(1);
		
		flip=rn_stream::throw_coin();
		size_type pa=flip?0:1;
		size_type pb=flip?1:0;
		
		value_type Ehat1=std::sqrt(P1*P1+this->m2(pa));
		value_type Ehat2=std::sqrt(P2*P2+this->m2(pb));

		this->p(pa).assign(0);
		this->p(pa)[k0]=Ehat1;
		this->p(pa)[kL]=P1;
		this->p(pb).assign(0);
		this->p(pb)[k0]=Ehat2;
		this->p(pb)[kL]=-P2;
		this->weight()=y_gen->weight();

		return true;
	    }

	    /* Weight evaluation implementation: */
	    
	    bool evaluate_weight()
	    {
		const momentum_type& p1=this->p(0);
		const momentum_type& p2=this->p(1);
		if(p1[k0]<this->m(0) or p2[k0]<this->m(1))
		{
		    this->weight()=(value_type)0;
                    this->integrand()=(value_type)0;
		    return false;
		}
		value_type P1=std::abs(p1[kL]);
		value_type P2=std::abs(p2[kL]);
		x1=P1/this->beam_energy(0);
		x2=P2/this->beam_energy(1);
		tau=x1*x2;

		value_type ybound=std::max(ymax,-(value_type)0.5*std::log(tau));
		y_gen->set_s_min(-ybound);
		y_gen->set_s_max(ybound);
		y=(value_type)0.5*std::log(x1/x2);
		if(!y_gen->evaluate_weight())
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		this->weight()=y_gen->weight();
		return true;
	    }

	    /* Flux factor output: */

	    value_type flux_factor() const
	    {
		size_type pa=flip?0:1;
		size_type pb=flip?1:0;
		int q1=(this->parton_id(pa)==21)?0:((hadron1<0)?(-this->parton_id(pa)):(this->parton_id(pa)));
		int q2=(this->parton_id(pb)==21)?0:((hadron2<0)?(-this->parton_id(pb)):(this->parton_id(pb)));
		
		return static_cast<value_type>(pdf_wrapper::ff(x1,x2,q1,q2,this->mu_F))*x1*x2/(this->s_hat()*std::sqrt(Kallen(this->s_hat(),this->m2(0),this->m2(1))));
	    }

	    /* Update method: */

	    void update()
	    {
		y_gen->integrand()=this->integrand();
		y_gen->update();
	    }
	    
	    /* Adaptor method implementation: */

	    void adapt_grids()
	    {
		y_gen->adapt();
	    }

	    /* Resets adaptive grids: */

	    void reset()
	    {
		this->base_type::reset();
		y_gen->reset();
	    }

	    /* Resets cross section: */

	    void reset_cross_section()
	    {
		this->base_type::reset_cross_section();
		y_gen->reset_cross_section();
	    }

	    /* Public readout methods */
	    /*------------------------*/

	    /* Clone method implementation: */

	    hadronic_is_y<model_t,rng_t,Minkowski_type>* clone() const
	    {
		return new hadronic_is_y<model_t,rng_t,Minkowski_type>(*this);
	    }

	    /* Beam particle id: */

	    int beam_id(size_type i) const
	    {
		if(i<2)
		{
		    return (i==0)?hadron1:hadron2;
		}
		return 0;
	    }

	    /* Cernlib pdf group number: */

	    int pdfg(size_type i) const
	    {
		return pdf_wrapper::group();
	    }

	    /* Cernlib pdf set number: */

	    int pdfs(size_type i) const
	    {
		return pdf_wrapper::set();
	    }

	    /* Returns the coupling constant corresponding to the factorization
	    scale: */

	    value_type alpha_s() const
	    {
		return pdf_wrapper::alpha_s(this->mu_F);
	    }

	    /* Serialisation: */
	    /*----------------*/

	    /* Returns the polymorphic type. */

	    std::string type() const
	    {
		switch(hadron1)
		{
		    case 2212:
			if(hadron2==2212)
			{
			    return "pp_y";
			}
			else if(hadron2==-2212)
			{
			    return "ppbar_y";
			}
			else
			{
			    return "X";
			}
			break;
		    case -2212:
			if(hadron2==2212)
			{
			    return "pbarp_y";
			}
			else if(hadron2==-2212)
			{
			    return "pbarpbar_y";
			}
			else
			{
			    return "X";
			}
			break;
		    default:
			return "X";
		}
	    }

	    /// Overloaded settings printing function.

	    std::ostream& print_settings(std::ostream& os) const
	    {
		os<<std::setw(30)<<std::left<<"beam energy -1:"<<this->beam_energy(0)<<std::endl;
		os<<std::setw(30)<<std::left<<"hadron -1:"<<hadron1<<std::endl;
		os<<std::setw(30)<<std::left<<"beam energy -2:"<<this->beam_energy(0)<<std::endl;
		os<<std::setw(30)<<std::left<<"hadron -2:"<<hadron2<<std::endl;
		os<<std::setw(30)<<std::left<<"Ecm:"<<this->Ecm()<<std::endl;
		os<<std::setw(30)<<std::left<<"PDF set:"<<pdf_name()<<std::endl;
		os<<std::setw(30)<<std::left<<"PDF nr:"<<pdf_number()<<std::endl;
		os<<std::setw(30)<<std::left<<"mu_F:"<<this->mu_F<<std::endl;
		os<<std::setw(30)<<std::left<<"alpha_s:"<<alpha_s()<<std::endl;
		return os;
	    }

	    /* Derived loading method: */

	    std::istream& load_data(std::istream& is)
	    {
		safe_read(is,tau_min);
		safe_read(is,ymax);
		safe_read(is,LHA_ymax);
		if(y_gen!=NULL)
		{
		    delete y_gen;
		}
		y_gen=s_generator<value_type,rng_t>::create_instance(&y,is,NULL,NULL,NULL);
		return is;
	    }

	    /* Derived saving method: */

	    std::ostream& save_data(std::ostream& os) const
	    {
		safe_write(os,tau_min);
		os<<"\t";
		safe_write(os,ymax);
		os<<"\t";
		safe_write(os,LHA_ymax);
		os<<std::endl;
		y_gen->save(os);
		return os;
	    }

	private:
	    
	    value_type y,tau,x1,x2,tau_min,LHA_ymax,ymax;
	    s_generator<value_type,rng_t>* y_gen;
	    bool flip;
    };
    template<class model_t,class rng_t>const typename hadronic_is_y<model_t,rng_t,Minkowski_type>::size_type hadronic_is_y<model_t,rng_t,Minkowski_type>::k0;
    template<class model_t,class rng_t>const typename hadronic_is_y<model_t,rng_t,Minkowski_type>::size_type hadronic_is_y<model_t,rng_t,Minkowski_type>::kL;
}

#endif /*CAMGEN_HAD_IS_H_*/

