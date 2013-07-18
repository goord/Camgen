//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PS_FAC_H_
#define CAMGEN_PS_FAC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * The phase space factory class connects particle types in Camgen with         *
 * propagator sampling types in Camgen and automates the creation of phase space *
 * branchings. The factories contain only static methods and data members.       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/ps_channel.h>
#include <Camgen/p_channel.h>
#include <Camgen/s_branch.h>
#include <Camgen/t_branch.h>
#include <Camgen/sum_branch.h>
#include <Camgen/s_branch_back.h>
#include <Camgen/t_branch_back.h>
#include <Camgen/sum_branch_back.h>

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t>class ps_factory
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_type::value_type value_type;
	    typedef spacetime_t spacetime_type;
	    typedef momentum_channel<model_t,N_in,N_out,rng_t> momentum_channel_type;
	    typedef particle_channel<model_t,N_in,N_out,rng_t> particle_channel_type;
	    typedef ps_branching<model_t,N_in,N_out,rng_t> branching_type;
	    typedef particle<model_t> particle_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;
	    typedef typename momentum_channel_type::momentum_type momentum_type;
	    typedef typename momentum_channel_type::size_type size_type;


	    /* Creates a particle channel instance from an input stream: */

	    static particle_channel_type* create_particle_channel(std::istream& is,momentum_channel_type* ps_channel_,const value_type* Ecmhat=NULL)
	    {
		std::string phi;
		value_type nu;
		bool r=false;
		is>>phi>>nu>>r;
		const particle_type* part=(phi=="X")?NULL:model_wrapper<model_t>::get_particle(phi);
		particle_channel_type* result=new particle_channel_type(ps_channel_,part,Ecmhat);
		result->nu=nu;
		result->redundant=r;
		delete (result->s_gen);
		bool adaptive=false;
		bool total_s=(ps_channel_->timelike() and ps_channel_->bitstring.count()==N_out);
		const value_type* m;
		if(total_s)
		{
		    if(Ecmhat!=NULL)
		    {
			m=Ecmhat;
		    }
		    else
		    {
			m=(part==NULL)?NULL:(part->get_mass_address());
		    }
		}
		else
		{
			m=(part==NULL)?NULL:(part->get_mass_address());
		}
		const value_type* w=(part==NULL)?NULL:(part->get_width_address());
		value_type* s=&(ps_channel_->invariant_mass);
		result->s_gen=s_generator<value_type,rng_t>::create_instance(s,is,m,w,&(result->nu));
		result->load(is);
		std::string endflag;
		while(endflag!="</particle>" and !is.eof())
		{
		    std::getline(is,endflag);
		}
		return result;
	    }

	    /* Factory method from input stream with allocated momentum: */

	    static momentum_channel_type* create_ps_channel(std::istream& is,const value_type* Ecmhat=NULL)
	    {
		bit_string_type bs;
		is>>bs;
		momentum_channel_type* result=new momentum_channel_type(bs);
		result->load(is);
		std::string flag;
		while(flag!="<particle>" and flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before particle channel data are read--returning NULL"<<endlog;
		    delete result;
		    return NULL;
		}
		particle_channel_type* phi;
		while(flag=="<particle>" and !is.eof())
		{
		    phi=create_particle_channel(is,result,Ecmhat);
		    if(phi!=NULL)
		    {
			result->particle_channels.push_back(phi);
		    }
		    std::getline(is,flag);
		}
		while(flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		return result;
	    }

	    /* Factory method from input stream with external momentum instance:
	     * */

	    static momentum_channel_type* create_ps_channel(momentum_type* p,std::istream& is,const value_type* Ecmhat=NULL)
	    {
		bit_string_type bs;
		is>>bs;
		momentum_channel_type* result=new momentum_channel_type(bs,p);
		result->load(is);
		std::string flag;
		while(flag!="<particle>" and flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read--returning NULL"<<endlog;
		    delete result;
		    return NULL;
		}
		particle_channel_type* phi;
		while(flag=="<particle>" and !is.eof())
		{
		    phi=create_particle_channel(is,result,Ecmhat);
		    if(phi!=NULL)
		    {
			result->particle_channels.push_back(phi);
		    }
		    std::getline(is,flag);
		}
		while(flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		return result;
	    }
	    

	    /* Creates a branching instance from the input stream: */

	    static branching_type* create_branching(std::istream& is,const vector<particle_channel_type*,N_in>& in_channels,const std::list<momentum_channel_type*>& ps_channels,const vector<particle_channel_type*,N_out>& out_channels)
	    {
		std::string br_id;
		is>>br_id;
		branching_type* result=NULL;
		if(br_id=="+" or br_id=="+~")
		{
		    std::string phi0,phi1,phi2;
		    is>>phi0>>phi1>>phi2;
		    particle_channel_type* c0=find_channel(phi0,in_channels,ps_channels,out_channels);
		    if(c0==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c1=find_channel(phi1,in_channels,ps_channels,out_channels);
		    if(c1==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c2=find_channel(phi2,in_channels,ps_channels,out_channels);
		    if(c2==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    if(br_id=="+")
		    {
			result=new sum_branching<model_t,N_in,N_out,rng_t>(c0,c1,c2);
		    }
		    else
		    {
			result=new backward_sum_branching<model_t,N_in,N_out,rng_t>(c0,c1,c2);
		    }
		    c0->branchings.push_back(result);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-generator id "<<br_id<<" not recognised--returning NULL."<<endlog;
		    move_to_end(is);
		    return NULL;
		}
		result->load(is);
		move_to_end(is);
		return result;
	    }

	private:

	    /* Static utility function for read-in from file: */

	    static particle_channel_type* find_channel(const std::string& name,const vector<particle_channel_type*,N_in>& in_channels,const std::list<momentum_channel_type*>& ps_channels,const vector<particle_channel_type*,N_out>& out_channels)
	    {
		for(size_type i=0;i<N_in;++i)
		{
		    if(name==in_channels[i]->name)
		    {
			return in_channels[i];
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    if(name==out_channels[i]->name)
		    {
			return out_channels[i];
		    }
		}
		for(typename std::list<momentum_channel_type*>::iterator it=ps_channels.begin();it!=ps_channels.end();++it)
		{
		    for(typename momentum_channel_type::particle_channel_list::const_iterator it2=(*it)->particle_channels.begin();it2!=(*it)->particle_channels.end();++it2)
		    {
			if(name==(*it2)->name)
			{
			    return *it2;
			}
		    }
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle channel "<<name<<" not found in list--returning NULL"<<endlog;
		return NULL;
	    }

	    /* Another static utility function for serialization: */

	    static std::istream& move_to_end(std::istream& is)
	    {
		std::string endflag;
		while(endflag!="</branching>" and !is.eof())
		{
		    std::getline(is,endflag);
		}
		return is;
	    }
    };

    template<class model_t,std::size_t N_out,class rng_t>class ps_factory<model_t,1,N_out,rng_t,Minkowski_type>
    {
	public:

	    /* Type definitions: */

	    typedef model_t model_type;
	    typedef typename model_type::value_type value_type;
	    typedef typename Minkowski_type::template implementation<typename model_t::value_type,model_t::dimension> spacetime_type;
	    typedef momentum_channel<model_t,1,N_out,rng_t> momentum_channel_type;
	    typedef particle_channel<model_t,1,N_out,rng_t> particle_channel_type;
	    typedef ps_branching<model_t,1,N_out,rng_t> branching_type;
	    typedef particle<model_t> particle_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;
	    typedef typename momentum_channel_type::momentum_type momentum_type;
	    typedef typename momentum_channel_type::size_type size_type;

	    /* Creates a forward s-type branching of the channel in. */

	    static branching_type* s_branch(particle_channel_type* in,particle_channel_type* out1,particle_channel_type* out2)
	    {
		bit_string_type bm=(out1->bitstring() & out2->bitstring());
		bit_string_type bp=(out1->bitstring() | out2->bitstring());
		if(bm.any() or bp!=in->bitstring())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-branching into invalid channels "<<out1->bitstring()<<","<<out2->bitstring()<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		return new s_branching<model_t,1,N_out,rng_t,spacetime_type>(in,out1,out2);
	    }

	    /* Creates a backward s-type branching of the channel in. */

	    static branching_type* backward_s_branch(particle_channel_type* in,particle_channel_type* out1,particle_channel_type* out2)
	    {
		bit_string_type bm=(out1->bitstring() & out2->bitstring());
		bit_string_type bp=(out1->bitstring() | out2->bitstring());
		if(bm.any() or bp!=in->bitstring())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-branching into invalid channels "<<out1->bitstring()<<","<<out2->bitstring()<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		return new backward_s_branching<model_t,1,N_out,rng_t,spacetime_type>(in,out1,out2);
	    }

	    /* Creates a particle channel instance from an input stream: */

	    static particle_channel_type* create_particle_channel(std::istream& is,momentum_channel_type* ps_channel_,const value_type* Ecmhat=NULL)
	    {
		std::string phi;
		value_type nu;
		bool r=false;
		is>>phi>>nu>>r;
		const particle_type* part=(phi=="X")?NULL:model_wrapper<model_t>::get_particle(phi);
		particle_channel_type* result=new particle_channel_type(ps_channel_,part,Ecmhat);
		result->nu=nu;
		result->redundant=r;
		delete (result->s_gen);
		bool total_s=(ps_channel_->bitstring.count()==N_out);
		const value_type* m;
		if(total_s)
		{
		    if(Ecmhat!=NULL)
		    {
			m=Ecmhat;
		    }
		    else
		    {
			m=(part==NULL)?NULL:(part->get_mass_address());
		    }
		}
		else
		{
			m=(part==NULL)?NULL:(part->get_mass_address());
		}
		const value_type* w=(part==NULL)?NULL:(part->get_width_address());
		value_type* s=&(ps_channel_->invariant_mass);
		result->s_gen=s_generator<value_type,rng_t>::create_instance(s,is,m,w,&(result->nu));
		result->load(is);
		std::string endflag;
		while(endflag!="</particle>" and !is.eof())
		{
		    std::getline(is,endflag);
		}
		return result;
	    }

	    /* Factory method from input stream with allocated momentum: */

	    static momentum_channel_type* create_ps_channel(std::istream& is,const value_type* Ecmhat=NULL)
	    {
		bit_string_type bs;
		is>>bs;
		momentum_channel_type* result=new momentum_channel_type(bs);
		result->load(is);
		std::string flag;
		while(flag!="<particle>" and flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before particle channel data are read--returning NULL"<<endlog;
		    delete result;
		    return NULL;
		}
		particle_channel_type* phi;
		while(flag=="<particle>" and !is.eof())
		{
		    phi=create_particle_channel(is,result,Ecmhat);
		    if(phi!=NULL)
		    {
			result->particle_channels.push_back(phi);
		    }
		    std::getline(is,flag);
		}
		while(flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		return result;
	    }

	    /* Factory method from input stream with external momentum instance:
	     * */

	    static momentum_channel_type* create_ps_channel(momentum_type* p,std::istream& is,const value_type* Ecmhat=NULL)
	    {
		bit_string_type bs;
		is>>bs;
		momentum_channel_type* result=new momentum_channel_type(bs,p);
		result->load(is);
		std::string flag;
		while(flag!="<particle>" and flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before initial data are read--returning NULL"<<endlog;
		    delete result;
		    return NULL;
		}
		particle_channel_type* phi;
		while(flag=="<particle>" and !is.eof())
		{
		    phi=create_particle_channel(is,result,Ecmhat);
		    if(phi!=NULL)
		    {
			result->particle_channels.push_back(phi);
		    }
		    std::getline(is,flag);
		}
		while(flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		return result;
	    }

	    /* Creates a branching instance from the input stream: */

	    static ps_branching<model_t,1,N_out,rng_t>* create_branching(std::istream& is,particle_channel_type* in_channel,const std::list<momentum_channel_type*>& ps_channels,const vector<particle_channel_type*,N_out>& out_channels)
	    {
		std::string br_id;
		is>>br_id;
		branching_type* result=NULL;
		if(br_id=="s" or br_id=="s~")
		{
		    std::string phi0,phi1,phi2;
		    is>>phi0>>phi1>>phi2;
		    particle_channel_type* c0=find_channel(phi0,in_channel,ps_channels,out_channels);
		    if(c0==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c1=find_channel(phi1,in_channel,ps_channels,out_channels);
		    if(c1==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c2=find_channel(phi2,in_channel,ps_channels,out_channels);
		    if(c2==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    if(br_id=="s")
		    {
			result=new s_branching<model_t,1,N_out,rng_t,spacetime_type>(c0,c1,c2);
		    }
		    else
		    {
			result=new backward_s_branching<model_t,1,N_out,rng_t,spacetime_type>(c0,c1,c2);
		    }
		    c0->branchings.push_back(result);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-generator id "<<br_id<<" not recognised--returning NULL"<<endlog;
		    move_to_end(is);
		    return NULL;
		}
		result->load(is);
		move_to_end(is);
		return result;
	    }

	private:

	    /* Static utility function for read-in from file: */

	    static particle_channel_type* find_channel(const std::string& name,particle_channel_type* in_channel,const std::list<momentum_channel_type*>& ps_channels,const vector<particle_channel_type*,N_out>& out_channels)
	    {
		if(name==in_channel->name)
		{
		    return in_channel;
		}
		for(size_type i=0;i<N_out;++i)
		{
		    if(name==out_channels[i]->name)
		    {
			return out_channels[i];
		    }
		}
		for(typename std::list<momentum_channel_type*>::const_iterator it=ps_channels.begin();it!=ps_channels.end();++it)
		{
		    for(typename momentum_channel_type::particle_channel_list::const_iterator it2=(*it)->particle_channels.begin();it2!=(*it)->particle_channels.end();++it2)
		    {
			if(name==(*it2)->name)
			{
			    return *it2;
			}
		    }
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle channel "<<name<<" not found in list--returning NULL"<<endlog;
		return NULL;
	    }

	    /* Another static utility function for serialization: */

	    static std::istream& move_to_end(std::istream& is)
	    {
		std::string endflag;
		while(endflag!="</branching>" and !is.eof())
		{
		    std::getline(is,endflag);
		}
		return is;
	    }
    };

    /* Phase space factory implementation for 2->N scattering processes. */

    template<class model_t,std::size_t N_out,class rng_t>class ps_factory<model_t,2,N_out,rng_t,Minkowski_type>
    {
	public:

	    typedef model_t model_type;
	    typedef typename model_type::value_type value_type;
	    typedef typename Minkowski_type::template implementation<typename model_t::value_type,model_t::dimension> spacetime_type;
	    typedef momentum_channel<model_t,2,N_out,rng_t> momentum_channel_type;
	    typedef particle_channel<model_t,2,N_out,rng_t> particle_channel_type;
	    typedef std::multimap< particle_channel_type*,particle_channel_type*> particle_channel_map;
	    typedef ps_branching<model_t,2,N_out,rng_t> branching_type;
	    typedef particle<model_t> particle_type;
	    typedef typename momentum_channel_type::bit_string_type bit_string_type;
	    typedef typename momentum_channel_type::momentum_type momentum_type;
	    typedef typename momentum_channel_type::size_type size_type;

	    /* Creates a forward s-type branching of the channel in: */

	    static branching_type* s_branch(particle_channel_type* in,particle_channel_type* out1,particle_channel_type* out2)
	    {
		if(in->bitstring()[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"spacelike incoming channel "<<in->bitstring()<<" allows no s-branching--returning NULL"<<endlog;
		    return NULL;
		}
		bit_string_type bm=(out1->bitstring() & out2->bitstring());
		bit_string_type bp=(out1->bitstring() | out2->bitstring());
		if(bm.any() or bp!=in->bitstring())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-branching into invalid channels "<<out1->bitstring()<<","<<out2->bitstring()<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		return new s_branching<model_t,2,N_out,rng_t,spacetime_type>(in,out1,out2);
	    }

	    /* Creates a backward s-type branching of the channel in: */

	    static branching_type* backward_s_branch(particle_channel_type* in,particle_channel_type* out1,particle_channel_type* out2)
	    {
		if(in->bitstring()[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"spacelike incoming channel "<<in->bitstring()<<" allows no s-branching--returning NULL"<<endlog;
		    return NULL;
		}
		bit_string_type bm=(out1->bitstring() & out2->bitstring());
		bit_string_type bp=(out1->bitstring() | out2->bitstring());
		if(bm.any() or bp!=in->bitstring())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-branching into invalid channels "<<out1->bitstring()<<","<<out2->bitstring()<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		return new backward_s_branching<model_t,2,N_out,rng_t,spacetime_type>(in,out1,out2);
	    }

	    /* Creates a forward t-type branching of the channel in1: */
	    
	    static branching_type* t_branch(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out1,particle_channel_type* out2,particle_channel_type* out3)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type s1=out1->bitstring();
		bit_string_type st=out2->bitstring();
		bit_string_type s2=out3->bitstring();

		if(!sa[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"incoming channel "<<sa<<" cannot be part of a t-channel--returning NULL"<<endlog;
		    return NULL;
		}
		if(!(sb[0] and sb.count()==1))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming momentum "<<sb<<" is not the second beam momentum--returning NULL"<<endlog;
		    return NULL;
		}
		if(s1[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first outgoing channel "<<s1<<" is not timelike--returning NULL"<<endlog;
		    return NULL;
		}
		if((st & s1).any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		if((st | s1)!=sa)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" and incoming channel "<<sa<<" ignored--returning NULL"<endlog;
		    return NULL;
		}
		bit_string_type bt=st;
		bt.reset(0);
		if(bt!=s2)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and second s-channel "<<s2<<" ignored--returning NULL"<<endlog;
		   return NULL;
		}
		return new t_branching<model_t,N_out,rng_t,spacetime_type>(in1,in2,out1,out2,out3,false);
	    }

	    /* Creates a t-type branching of the channel in1: */
	    
	    static branching_type* backward_t_branch(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out1,particle_channel_type* out2,particle_channel_type* out3)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type s1=out1->bitstring();
		bit_string_type st=out2->bitstring();
		bit_string_type s2=out3->bitstring();

		if(!sa[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"incoming channel "<<sa<<" cannot be part of a t-channel--returning NULL"<<endlog;
		    return NULL;
		}
		if(!(sb[0] and sb.count()==1))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming momentum "<<sb<<" is not the second beam momentum--returning NULL"<<endlog;
		    return NULL;
		}
		if(s1[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first outgoing channel "<<s1<<" is not timelike--returning NULL"<<endlog;
		    return NULL;
		}
		if((st & s1).any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" ignored--returning NULL"<<endlog;
		    return NULL;
		}
		if((st | s1)!=sa)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" and incoming channel "<<sa<<" ignored--returning NULL"<endlog;
		    return NULL;
		}
		bit_string_type bt=st;
		bt.reset(0);
		if(bt!=s2)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and second s-channel "<<s2<<" ignored--returning NULL"<<endlog;
		   return NULL;
		}
		return new backward_t_branching<model_t,N_out,rng_t,spacetime_type>(in1,in2,out1,out2,out3,false);
	    }

	    /* Creates a forward momentum sum branching: */

	    static branching_type* sum_branch(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type ss=out->bitstring();
		sa.flip();
		
		if(sa.any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first incoming channel "<<sa<<" is not the first beam--returning NULL"<<endlog;
		    return NULL;
		}
		if(!(sb.count()==1 and sb[0]))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming channel "<<sb<<" is not the second beam--returning NULL"<<endlog;
		    return NULL;
		}
		if(ss[0] or ss.count()!=N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"outgoing channel "<<ss<<" is not the total outgoing momentum--returning NULL"<<endlog;
		    return NULL;
		}
		return new sum_branching<model_t,2,N_out,rng_t>(in1,in2,out);
	    }

	    /* Creates a backward momentum sum branching: */

	    static branching_type* backward_sum_branch(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type ss=out->bitstring();
		sa.flip();
		
		if(sa.any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first incoming channel "<<sa<<" is not the first beam--returning NULL"<<endlog;
		    return NULL;
		}
		if(!(sb.count()==1 and sb[0]))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming channel "<<sb<<" is not the second beam--returning NULL"<<endlog;
		    return NULL;
		}
		if(ss[0] or ss.count()!=N_out)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"outgoing channel "<<ss<<" is not the total outgoing momentum--returning NULL"<<endlog;
		    return NULL;
		}
		return new backward_sum_branching<model_t,2,N_out,rng_t>(in1,in2,out);
	    }

	    /* Registers the first channel as a t-channel and the second channel
	     * as the second outgoing s-channel in the corresponding 2->2
	     * t-diagram: */

	    static bool register_t_channel(particle_channel_type* in,particle_channel_type* out,bool final)
	    {
		bit_string_type sa=in->bitstring();
		bit_string_type ss=out->bitstring();
		if(!sa[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first incoming channel "<<sa<<" cannot be part of a t-channel"<<endlog;
		    return false;
		}
		bit_string_type bt=sa;
		bt.reset(0);
		if(bt!=ss)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<sa<<" and second s-channel "<<ss<<" ignored"<<endlog;
		   return false;
		}
		if(in->bitstring().count()<=N_out)
		{
		    if(final)
		    {
			final_t_channels.insert(std::pair<particle_channel_type*,particle_channel_type*>(in,out));
		    }
		    else
		    {
			inter_t_channels.insert(std::pair<particle_channel_type*,particle_channel_type*>(in,out));
		    }
		    return true;
		}
		return false;
	    }

	    /* Creates all forward t-type branchings, looking up the second s-channel in
	    the register entries corresponding to the t-channel out2. */

	    static std::vector<branching_type*> t_branch(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out1,particle_channel_type* out2)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type s1=out1->bitstring();
		bit_string_type st=out2->bitstring();

		std::vector<branching_type*>result;

		if(!sa[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"incoming channel "<<sa<<" cannot be part of a t-channel"<<endlog;
		    return result;
		}
		if(!(sb[0] and sb.count()==1))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming momentum "<<sb<<" is not the second beam momentum"<<endlog;
		    return result;
		}
		if(s1[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first outgoing channel "<<s1<<" is not timelike"<<endlog;
		    return result;
		}
		if((st & s1).any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" ignored"<<endlog;
		    return result;
		}
		if((st | s1)!=sa)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" and incoming channel "<<sa<<" ignored"<<endlog;
		    return result;
		}
		bit_string_type bt=st;
		bt.reset(0);
		for(typename particle_channel_map::iterator it=final_t_channels.equal_range(out2).first;it!=final_t_channels.equal_range(out2).second;++it)
		{
		    if(bt!=it->second->bitstring())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and second s-channel "<<(it->second->bitstring())<<" ignored"<<endlog;
			continue;
		    }
		    branching_type* br=new t_branching<model_t,N_out,rng_t,spacetime_type>(in1,in2,out1,out2,it->second,true);
		    result.push_back(br);
		}
		for(typename particle_channel_map::iterator it=inter_t_channels.equal_range(out2).first;it!=inter_t_channels.equal_range(out2).second;++it)
		{
		    if(bt!=it->second->bitstring())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and second s-channel "<<(it->second->bitstring())<<" ignored"<<endlog;
			continue;
		    }
		    branching_type* br=new t_branching<model_t,N_out,rng_t,spacetime_type>(in1,in2,out1,out2,it->second,false);
		    result.push_back(br);
		}
		return result;
	    }

	    /* Creates all t-type branchings, looking up the second s-channel in
	    the register entries corresponding to the t-channel out2. */

	    static std::vector<branching_type*> backward_t_branch(particle_channel_type* in1,particle_channel_type* in2,particle_channel_type* out1,particle_channel_type* out2)
	    {
		bit_string_type sa=in1->bitstring();
		bit_string_type sb=in2->bitstring();
		bit_string_type s1=out1->bitstring();
		bit_string_type st=out2->bitstring();

		std::vector<branching_type*>result;

		if(!sa[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"incoming channel "<<sa<<" cannot be part of a t-channel"<<endlog;
		    return result;
		}
		if(!(sb[0] and sb.count()==1))
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"second incoming momentum "<<sb<<" is not the second beam momentum"<<endlog;
		    return result;
		}
		if(s1[0])
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first outgoing channel "<<s1<<" is not timelike"<<endlog;
		    return result;
		}
		if((st & s1).any())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" ignored"<<endlog;
		    return result;
		}
		if((st | s1)!=sa)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and first s-channel "<<s1<<" and incoming channel "<<sa<<" ignored"<<endlog;
		    return result;
		}
		bit_string_type bt=st;
		bt.reset(0);
		for(typename particle_channel_map::iterator it=final_t_channels.equal_range(out2).first;it!=final_t_channels.equal_range(out2).second;++it)
		{
		    if(bt!=it->second->bitstring())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and second s-channel "<<(it->second->bitstring())<<" ignored"<<endlog;
			continue;
		    }
		    branching_type* br=new backward_t_branching<model_t,N_out,rng_t,spacetime_type>(in1,in2,out1,out2,it->second,true);
		    result.push_back(br);
		}
		for(typename particle_channel_map::iterator it=inter_t_channels.equal_range(out2).first;it!=inter_t_channels.equal_range(out2).second;++it)
		{
		    if(bt!=it->second->bitstring())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"invalid combination of t-channel "<<st<<" and second s-channel "<<(it->second->bitstring())<<" ignored"<<endlog;
			continue;
		    }
		    branching_type* br=new backward_t_branching<model_t,N_out,rng_t,spacetime_type>(in1,in2,out1,out2,it->second,false);
		    result.push_back(br);
		}
		return result;
	    }

	    /* Creates a particle channel instance from an input stream: */

	    static particle_channel_type* create_particle_channel(std::istream& is,momentum_channel_type* ps_channel_,const value_type* Ecmhat=NULL)
	    {
		std::string phi;
		value_type nu;
		bool r=false;
		is>>phi>>nu>>r;
		const particle_type* part=(phi=="X")?NULL:model_wrapper<model_t>::get_particle(phi);
		particle_channel_type* result=new particle_channel_type(ps_channel_,part,Ecmhat);
		result->nu=nu;
		result->redundant=r;
		delete (result->s_gen);
		bool total_s=(ps_channel_->timelike() and ps_channel_->bitstring.count()==N_out);
		const value_type* m;
		if(total_s)
		{
		    if(Ecmhat!=NULL)
		    {
			m=Ecmhat;
		    }
		    else
		    {
			m=(part==NULL)?NULL:(part->get_mass_address());
		    }
		}
		else
		{
			m=(part==NULL)?NULL:(part->get_mass_address());
		}
		const value_type* w=(part==NULL)?NULL:(part->get_width_address());
		value_type* s=&(ps_channel_->invariant_mass);
		result->s_gen=s_generator<value_type,rng_t>::create_instance(s,is,m,w,&(result->nu));
		result->load(is);
		std::string endflag;
		while(endflag!="</particle>" and !is.eof())
		{
		    std::getline(is,endflag);
		}
		return result;
	    }

	    /* Factory method from input stream with allocated momentum: */

	    static momentum_channel_type* create_ps_channel(std::istream& is,const value_type* Ecmhat=NULL)
	    {
		bit_string_type bs;
		is>>bs;
		momentum_channel_type* result=new momentum_channel_type(bs);
		result->load(is);
		std::string flag;
		while(flag!="<particle>" and flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before particle channel data are read--returning NULL"<<endlog;
		    delete result;
		    return NULL;
		}
		particle_channel_type* phi;
		while(flag=="<particle>" and !is.eof())
		{
		    phi=create_particle_channel(is,result,Ecmhat);
		    if(phi!=NULL)
		    {
			result->particle_channels.push_back(phi);
		    }
		    std::getline(is,flag);
		}
		while(flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		return result;
	    }

	    /* Factory method from input stream with external momentum instance:
	     * */

	    static momentum_channel_type* create_ps_channel(momentum_type* p,std::istream& is,const value_type* Ecmhat=NULL)
	    {
		bit_string_type bs;
		is>>bs;
		momentum_channel_type* result=new momentum_channel_type(bs,p);
		result->load(is);
		std::string flag;
		while(flag!="<particle>" and flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		if(is.eof())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"end of file reached before particle channel data are read--NULL returned"<<endlog;
		    delete result;
		    return NULL;
		}
		particle_channel_type* phi;
		while(flag=="<particle>" and !is.eof())
		{
		    phi=create_particle_channel(is,result,Ecmhat);
		    if(phi!=NULL)
		    {
			result->particle_channels.push_back(phi);
		    }
		    std::getline(is,flag);
		}
		while(flag!="</channel>" and !is.eof())
		{
		    std::getline(is,flag);
		}
		return result;
	    }

	    /* Creates a branching instance from the input stream: */

	    static branching_type* create_branching(std::istream& is,const vector<particle_channel_type*,2>& in_channels,const std::list<momentum_channel_type*>& ps_channels,const vector<particle_channel_type*,N_out>& out_channels)
	    {
		std::string br_id;
		is>>br_id;
		branching_type* result=NULL;
		if(br_id=="s" or br_id=="s~")
		{
		    std::string phi0,phi1,phi2;
		    is>>phi0>>phi1>>phi2;
		    particle_channel_type* c0=find_channel(phi0,in_channels,ps_channels,out_channels);
		    if(c0==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c1=find_channel(phi1,in_channels,ps_channels,out_channels);
		    if(c1==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c2=find_channel(phi2,in_channels,ps_channels,out_channels);
		    if(c2==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    if(br_id=="s")
		    {
			result=new s_branching<model_t,2,N_out,rng_t,spacetime_type>(c0,c1,c2);
		    }
		    else if(br_id=="s~")
		    {
			result=new backward_s_branching<model_t,2,N_out,rng_t,spacetime_type>(c0,c1,c2);
		    }
		    c0->branchings.push_back(result);
		}
		else if(br_id=="t" or br_id=="t." or br_id=="t~" or br_id=="t~.")
		{
		    std::string phi0,phi1,phi2,phi3,phi4;
		    is>>phi0>>phi1>>phi2>>phi3>>phi4;
		    particle_channel_type* c0=find_channel(phi0,in_channels,ps_channels,out_channels);
		    if(c0==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c1=find_channel(phi1,in_channels,ps_channels,out_channels);
		    if(c1==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c2=find_channel(phi2,in_channels,ps_channels,out_channels);
		    if(c2==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c3=find_channel(phi3,in_channels,ps_channels,out_channels);
		    if(c3==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c4=find_channel(phi4,in_channels,ps_channels,out_channels);
		    if(c4==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    if(br_id=="t" or br_id=="t.")
		    {
			bool final=(br_id=="t.");
			result=new t_branching<model_t,N_out,rng_t,spacetime_type>(c0,c1,c2,c3,c4,final);
		    }
		    else
		    {
			bool final=(br_id=="t~.");
			result=new backward_t_branching<model_t,N_out,rng_t,spacetime_type>(c0,c1,c2,c3,c4,final);
			std::string phi5;
			is>>phi5;
			particle_channel_type* c5=find_channel(phi5,in_channels,ps_channels,out_channels);
			static_cast<backward_t_branching<model_t,N_out,rng_t,spacetime_type>*>(result)->set_s12_channel(c5);
		    }
		    c0->branchings.push_back(result);
		}
		else if(br_id=="+" or br_id=="+~")
		{
		    std::string phi0,phi1,phi2;
		    is>>phi0>>phi1>>phi2;
		    particle_channel_type* c0=find_channel(phi0,in_channels,ps_channels,out_channels);
		    if(c0==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c1=find_channel(phi1,in_channels,ps_channels,out_channels);
		    if(c1==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    particle_channel_type* c2=find_channel(phi2,in_channels,ps_channels,out_channels);
		    if(c2==NULL)
		    {
			move_to_end(is);
			return NULL;
		    }
		    if(br_id=="+")
		    {
			result=new sum_branching<model_t,2,N_out,rng_t>(c0,c1,c2);
		    }
		    else
		    {
			result=new backward_sum_branching<model_t,2,N_out,rng_t>(c0,c1,c2);
		    }
		    c0->branchings.push_back(result);
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"s-generator id "<<br_id<<" not recognised--returning NULL"<<endlog;
		    move_to_end(is);
		    return NULL;
		}
		result->load(is);
		move_to_end(is);
		return result;
	    }

	    /* Clears the register of t-branchings: */

	    static void flush_t_channels()
	    {
		final_t_channels.clear();
		inter_t_channels.clear();
	    }

	private:

	    /* Utility map registering final t-channels. */

	    static particle_channel_map final_t_channels;

	    /* Utility map registering intermediate t-channels. */

	    static particle_channel_map inter_t_channels;

	    /* Static utility function for read-in from file: */

	    static particle_channel_type* find_channel(const std::string& name,const vector<particle_channel_type*,2>& in_channels,const std::list<momentum_channel_type*>& ps_channels,const vector<particle_channel_type*,N_out>& out_channels)
	    {
		for(size_type i=0;i<2;++i)
		{
		    if(name==in_channels[i]->name)
		    {
			return in_channels[i];
		    }
		}
		for(size_type i=0;i<N_out;++i)
		{
		    if(name==out_channels[i]->name)
		    {
			return out_channels[i];
		    }
		}
		for(typename std::list<momentum_channel_type*>::const_iterator it=ps_channels.begin();it!=ps_channels.end();++it)
		{
		    for(typename momentum_channel_type::particle_channel_list::iterator it2=(*it)->particle_channels.begin();it2!=(*it)->particle_channels.end();++it2)
		    {
			if(name==(*it2)->name)
			{
			    return *it2;
			}
		    }
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"particle channel "<<name<<" not found in list--returning NULL"<<endlog;
		return NULL;
	    }

	    /* Another static utility function for serialization: */

	    static std::istream& move_to_end(std::istream& is)
	    {
		std::string endflag;
		while(endflag!="</branching>" and !is.eof())
		{
		    std::getline(is,endflag);
		}
		return is;
	    }
    };
    template<class model_t,std::size_t N_out,class rng_t>typename ps_factory<model_t,2,N_out,rng_t,Minkowski_type>::particle_channel_map ps_factory<model_t,2,N_out,rng_t,Minkowski_type>::final_t_channels;
    template<class model_t,std::size_t N_out,class rng_t>typename ps_factory<model_t,2,N_out,rng_t,Minkowski_type>::particle_channel_map ps_factory<model_t,2,N_out,rng_t,Minkowski_type>::inter_t_channels;
}

#endif /*CAMGEN_PS_FAC_H_*/

