//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file input_ps.h
    \brief Specialisation for the phase space generator factory for the istream classes.
 */

#ifndef INPUT_PS_H_
#define INPUT_PS_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Phase space generator factory specialisation for the istream class, creating  *
 * an phase space generator instance from a filestream.                          *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/MC_config.h>
#include <Camgen/isgen_fac.h>
#include <Camgen/ps_tree.h>
#include <Camgen/rambo.h>

namespace Camgen
{
    /// Phase space generator factory.

    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_generator_factory
    {
	public:

	    /* Type definitions. */

	    typedef ps_generator<model_t,N_in,N_out> generator_type;
	    typedef phase_space_generators::type fs_generator_tag;
	    typedef initial_states::type is_generator_tag;

	    /// Factory method with recursive tree iterator.

	    static generator_type* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		return create_instance(it,initial_state_type(),phase_space_generator_type());
	    }

	    /// Factory method with recursive tree iterator, initial state
	    /// generator tag and phase space generator tag.

	    static generator_type* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it,is_generator_tag is_tag,fs_generator_tag fs_tag)
	    {
		bool backward_shat=(fs_tag==phase_space_generators::recursive_backward_shat);
		typename initial_state_factory<model_t,N_in,rng_t>::generator_type* isgen=initial_state_factory<model_t,N_in,rng_t>::create_instance(is_tag,backward_shat);
		if(isgen==NULL)
		{
		    return NULL;
		}
		fs_generator_tag fst=fs_tag;
		if(backward_shat and isgen->s_hat_sampling)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"initial state does not allow s-hat sampling by final state generator--continuing with IS s-hat sampling"<<endlog;
		    fst=phase_space_generators::recursive_backward_s;
		}
		generator_type* result=NULL;
		switch(fst)
		{
		    case phase_space_generators::uniform:
			result=rambo<model_t,N_in,N_out,rng_t>::create_instance(it,isgen);
			break;
		    case phase_space_generators::recursive:
			result=ps_tree<model_t,N_in,N_out,rng_t>::create_instance(it,isgen);
			break;
		    case phase_space_generators::recursive_backward_s:
			result=ps_tree<model_t,N_in,N_out,rng_t>::create_instance(it,isgen);
			break;
		    case phase_space_generators::recursive_backward_shat:
			result=ps_tree<model_t,N_in,N_out,rng_t>::create_instance(it,isgen);
			break;
		    default:
			break;
		}
		if(result!=NULL)
		{
		    for(int i=1;i<=(int)N_in;++i)
		    {
			result->set_beam_energy(-i,beam_energy(i));
		    }
		    result->refresh_Ecm();
		    for(int i=1;i<=(int)N_out;++i)
		    {
			for(int j=i+1;j<=(int)N_out;++j)
			{
			    result->set_m_min(i,j,basic_cuts::m_min(i,j));
			}
		    }
		    result->refresh_m_min();
		}
		return result;
	    }

	    /// Factory method with recursive tree iterator and input stream.
	    
	    static generator_type* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it,std::istream& is)
	    {
		std::string initline,endline;
		do
		{
		    std::getline(is,initline);
		}
		while(initline!="<psgen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read--returning NULL"<<endlog;
		    return NULL;
		}
		typename initial_state_factory<model_t,N_in,rng_t>::generator_type* isgen=initial_state_factory<model_t,N_in,rng_t>::create_instance(is);
		if(isgen==NULL)
		{
		    do
		    {
			std::getline(is,endline);
		    }
		    while(endline!="</psgen>" and !is.eof());
		    return NULL;
		}
		do
		{
		    std::getline(is,initline);
		}
		while(initline!="<fsgen>" and !is.eof());
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read--returning NULL"<<endlog;
		    return NULL;
		}
		generator_type* result;
		std::string fstype;
		is>>fstype;
		if(fstype=="rambo")
		{
		    result=new rambo<model_t,N_in,N_out,rng_t>(isgen);
		}
		else if(fstype=="pstree")
		{
		    bool q1=false,q2=false,q3=false;
		    typename generator_type::size_type b;
		    is>>q1>>q2>>q3>>b;
		    result=new ps_tree<model_t,N_in,N_out,rng_t>(isgen,false,q1,q2,q3,b);
		}
		else if(fstype=="pstree~")
		{
		    bool q1=false,q2=false,q3=false;
		    typename generator_type::size_type b;
		    is>>q1>>q2>>q3>>b;
		    result=new ps_tree<model_t,N_in,N_out,rng_t>(isgen,true,q1,q2,q3,b);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"momentum generator type "<<fstype<<" not recognised for "<<N_in<<" -> "<<N_out<<" process"<<endlog;
		    result=NULL;
		}
		if(result!=NULL)
		{
		    if(!result->set_tree(it))
		    {
			return NULL;
		    }
		    result->load(is);
		    result->refresh_Ecm();
		    result->refresh_m_min();
		}
		do
		{
		    std::getline(is,endline);
		}
		while(endline!="</psgen>" and !is.eof());
		return result;
	    }
    };
}

#endif /*INPUT_PS_H_*/

