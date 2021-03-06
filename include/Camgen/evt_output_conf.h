//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file evt_output_conf.h
  \brief Abstract base class for event output configurations.
  */

#ifndef CAMGEN_EVT_OUTPUT_CONF_H_
#define CAMGEN_EVT_OUTPUT_CONF_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract for base class for event output configurations, which should implement *
 * the fill() method to evaluate the variables that need to be written to disk.    *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <sstream>
#include <Camgen/evt_output.h>

namespace Camgen
{
    /// Event output configuration class. Implementors should add variables (branches in the output) and fill the
    /// corresponding pointers with event kinematic variables.

    template<class model_t,std::size_t N_in,std::size_t N_out>class event_output_configuration
    {
	public:

	    /* Type definitions: */

	    typedef event<model_t,N_in,N_out> event_type;
	    typedef typename event_type::momentum_type momentum_type;
	    typedef typename event_type::value_type value_type;
	    typedef typename event_type::size_type size_type;

	    /* Public constructors: */
	    /*----------------------*/

	    /// Default constructor.

	    event_output_configuration():output(NULL),evt_size(0){}

	    /// Constructor with generator and output instance arguments.

	    event_output_configuration(event_output<model_t,N_in,N_out>* output_):output(output_),evt_size(0){}

	    /* Public destructors: */
	    /*---------------------*/

	    virtual ~event_output_configuration(){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Copies the output variable addresses to the output object

	    void add_branches(event_output<model_t,N_in,N_out>* output_)
	    {
		output=output_;
		add_variables();
	    }

	    /// Abstract clone method.

	    virtual event_output_configuration<model_t,N_in,N_out>* clone() const=0;

	    /// Abstract method filling the output variables.

	    virtual void fill(const event_type&)=0;

	    /* Public readout functions: */
	    /*---------------------------*/

	    size_type event_size() const
	    {
		return evt_size;
	    }

	protected:

	    /// Abstract mathod adding variables to the interface. Derived classes should
	    /// implement this function with add_variable() statements to stream the
	    /// variable addresses to the output interface.

	    virtual void add_variables()=0;

	    /// Adds a momentum-valued output variable branch to the output tree.

	    bool add_variable(const momentum_type& q,const std::string& name)
	    {
		if(output!=NULL)
		{
		    evt_size+=sizeof(momentum_type);
		    return output->branch(&q,name);
		}
		return false;
	    }

	    /// Adds a floating-point output variable branch to the output tree.

	    bool add_variable(const value_type& q,const std::string& name)
	    {
		if(output!=NULL)
		{
		    evt_size+=sizeof(value_type);
		    return output->branch(&q,name);
		}
		return false;
	    }

	    /// Adds an integer-valued output variable branch to the output tree.

	    bool add_variable(const int& q,const std::string& name)
	    {
		if(output!=NULL)
		{
		    evt_size+=sizeof(int);
		    return output->branch(&q,name);
		}
		return false;
	    }

	    /// Adds an boolean-valued output variable branch to the output tree.

	    bool add_variable(const bool& q,const std::string& name)
	    {
		if(output!=NULL)
		{
		    evt_size+=sizeof(bool);
		    return output->branch(&q,name);
		}
		return false;
	    }

	private:

	    /* Output file object: */

	    event_output<model_t,N_in,N_out>* output;

	    /* event size: */

	    size_type evt_size;
    };

    /// Standard output configuration class. Writes all momenta plus partonic CM energy to the output stream.

    template<class model_t,std::size_t N_in,std::size_t N_out>class standard_output_configuration: public event_output_configuration<model_t,N_in,N_out>
    {
        typedef event_output_configuration<model_t,N_in,N_out> base_type;

        public:

	    /* Type definitions: */

	    typedef typename base_type::event_type event_type;
	    typedef typename base_type::momentum_type momentum_type;
	    typedef typename base_type::value_type value_type;
	    typedef typename base_type::size_type size_type;

	    /// Clone implementation.

	    event_output_configuration<model_t,N_in,N_out>* clone() const
            {
                return new standard_output_configuration<model_t,N_in,N_out>();
            }

            /// Implementation of the fill function. Copies Ecmhat and all momenta.

            void fill(const event_type& evt)
            {
                ecmhat=evt.Ecm_hat();
                for(size_type i=0;i<N_in;++i)
                {
                    pin[i]=evt.p_in(i);
                }
                for(size_type i=0;i<N_out;++i)
                {
                    pin[i]=evt.p_out(i);
                }
            }

        protected:

            /// Variable-adding function. Adds Ecmhat and all omenta to the output tree.

            void add_variables()
            {
                this->add_variable(ecmhat,"Ecmhat");
                for(size_type i=0;i<N_in;++i)
                {
                    std::stringstream ss;
                    ss<<"p"<<(i+1);
                    this->add_variable(pin[i],ss.str());
                }
                for(size_type i=0;i<N_out;++i)
                {
                    std::stringstream ss;
                    ss<<"q"<<(i+1);
                    this->add_variable(pout[i],ss.str());
                }
            }

        private:

            /* Local copy of the partonic invariant mass: */

            value_type ecmhat;

            /* Local copy of the incoming momenta: */

            vector<momentum_type,N_in> pin;
            
            /* Local copy of the outgoing momenta: */
            
            vector<momentum_type,N_out> pout;
    };
}

#endif /*CAMGEN_EVT_OUTPUT_CONF_H_*/

