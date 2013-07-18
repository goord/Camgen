//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_CURRENT_TREE_H_
#define CAMGEN_CURRENT_TREE_H_

#include <algorithm>
#include <Camgen/current.h>
#include <Camgen/logstream.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of the current_tree class, which is essentially a static list of   *
 * current objects. Given the model type, the number of in- and outgoing legs of *
 * the process and the final particle number, the initialisation of the data     *
 * tree creates a big vector of all possible currents (all particle flavours in  *
 * all possible momentum channels). The actual subprocess trees consist of lists *
 * of interactions between these currents.                                       *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out>class current_tree
    {
	public:

	    /* Declaring the interaction tree a friend: */

	    friend class process_tree<model_t,N_in,N_out>;

	    /* Defining the number of particles the trees are built upon: */

	    static const std::size_t N_bits=N_in+N_out-1;
	    
	    /* Number of momentum channels: */
	    
	    static const std::size_t static_size=(1<<N_bits);

	private:

	    /* Integer denoting the final particle for process trees built upon
	     * the current tree: */

	    static std::size_t N_final;

	public:

	    /* Current type definition: */

	    typedef current<model_t,N_bits,get_colour_treatment<model_t>::decomposes> current_type;
	    
	    /* Current iterator type definition: */
	    
	    typedef typename std::vector<current_type>::iterator iterator;
	    
	    /* Current vector size type definition: */
	    
	    typedef typename std::vector<current_type>::size_type size_type;

	    /* Initialisation phase: */

	    static void initialise(std::size_t N_f=0)
	    {
		if(!initialised and N_f<(N_in+N_out))
		{
		    model<model_t>::get_instance();
		    N_final=N_f;
		    refresh();
		    initialised=true;
		}
		else if(N_f!=N_final and N_f<(N_in+N_out))
		{
		    set_final_current(N_f);
		}
	    }

	    /* Tree construction function: */

	    static void refresh()
	    {
		if(initialised)
		{
		    log(log_level::message)<<"rebuilding current tree..."<<endlog;
		    log(log_level::warning)<<"previously built process trees will become invalid..."<<endlog;
		}

		/* Resize the current vector to N_channels x N_flavours: */

		flavours=model_wrapper<model_t>::flavours();
		data.clear();
		data.resize(static_size*(flavours));
		
		/* Construct initial & internal currents: */
		
		bit_string<N_bits>B;
		B.set(0);
		size_type counter=0;
		size_type m=(N_final<N_in)?(N_in-1):N_in;
		for(size_type i=0;i<static_size-1;++i)
		{
		    for(size_type j=0;j<flavours;++j)
		    {
			data[counter]=current_type(model_wrapper<model_t>::get_particle(j),B,(i>=m and i<N_bits));
			++counter;
		    }
		    B.next();
		}

		/* Construct final currents: */

		for(size_type j=0;j<flavours;++j)
		{
		    data[counter]=current_type(model_wrapper<model_t>::get_particle(j),B,(N_final>=N_in));
		    ++counter;
		}
	    }

	    /* Function re-assigning the final particle currents: */

	    static void set_final_current(size_type n)
	    {
		bool p=(N_final<N_in);
		bool q=(n<N_in);
		if(p and !q)
		{
		    log(log_level::message)<<"re-assigning momentum directions in current tree..."<<endlog;
		    log(log_level::warning)<<"previously built process trees will become invalid..."<<endlog;
		    
		    for(size_type i=(N_in-1)*flavours;i<N_in*flavours;++i)
		    {
			data[i].set_incoming();
		    }
		    for(size_type i=data.size()-flavours;i<data.size();++i)
		    {
			data[i].set_outgoing();
		    }
		}
		if(!p and q)
		{
		    log(log_level::message)<<"re-assigning momentum directions in current tree..."<<endlog;
		    log(log_level::warning)<<"previously built process trees will become invalid..."<<endlog;

		    for(size_type i=(N_in-1)*flavours;i<N_in*flavours;++i)
		    {
			data[i].set_outgoing();
		    }
		    for(size_type i=data.size()-flavours;i<data.size();++i)
		    {
			data[i].set_incoming();
		    }
		}
		N_final=n;
	    }

	    /* Final current block readout: */

	    static size_type final_current()
	    {
		return N_final;
	    }

	    /* Output method: */

	    static std::ostream& print(std::ostream& os)
	    {
		for(int i=0;i<data.size();++i)
		{
		    data[i].print(os);
		    os<<std::endl;
		}
		return os;
	    }

	    /* Useful iterators in the tree: */

	    static iterator begin()
	    {
		return data.begin();
	    }
	    static iterator end()
	    {
		return data.end();
	    }
	    static iterator begin_internals()
	    {
		if(initialised)
		{
		    return data.begin()+N_bits*flavours;
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"internal current iterator requested for uninitialised data"<<endlog;
		    return data.end();
		}
	    }
	    static iterator end_internals()
	    {
		if(initialised)
		{
		    return data.end()-flavours;
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"internal current iterator requested for uninitialised data"<<endlog;
		return data.end();
	    }

	    /* Current-finding algorithms: */

	    static iterator find_current(const bit_string<N_bits>& B)
	    {
		if(initialised)
		{
		    return data.begin()+((B.to_integer()-1)*flavours);
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"internal current iterator requested for uninitialised data"<<endlog;
		return data.end();
	    }

	    static iterator find_current(const bit_string<N_bits>& B,const particle<model_t>* phi)
	    {
		if(initialised)
		{
		    CAMGEN_ERROR_IF((phi==NULL),"attempt to dereference NULL particle type instance...");
		    return data.begin()+((B.to_integer()-1)*flavours+phi->get_flavour());
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"current iterator requested for uninitialised data"<<endlog;
		return data.end();
	    }

	    static iterator find_current(const bit_string<N_bits>& B,const std::string& str)
	    {
		return find_current(B,model_wrapper<model_t>::get_particle(str));
	    }

	    static iterator find_current(const bit_string<N_bits>& B,const std::size_t flav)
	    {
		if(initialised)
		{
		    return (flav<flavours)?(data.begin()+(B.to_integer()-1)*flavours+flav):data.end();
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"current iterator requested for uninitialised data"<<endlog;
		return data.end();
	    }

	    /* Finding the first marked current in momentum channel B: */

	    static iterator first_marked_current(const bit_string<N_bits>& B)
	    {
		if(initialised)
		{
		    iterator it=data.begin()+((B.to_integer()-1)*flavours);
		    for(size_type i=0;i<flavours;++i)
		    {
			if(it->is_marked())
			{
			    return it;
			}
			++it;
		    }
		}
		else
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"first marked current iterator requested for uninitialised data"<<endlog;
		}
		return data.end();
	    }

	    /* Finding the next marked current in the momentum channel of the
	     * argument iterator, and moving it to this place. The return value
	     * denotes whether the iterator shifted from the end of the channel
	     * block to the beginning. */

	    static bool next_marked_current(iterator& it)
	    {
		if(initialised)
		{
		    if(it>=data.end())
		    {
			log(log_level::warning)<<CAMGEN_STREAMLOC<<"argument iterator out of range"<<endlog;
			return false;
		    }
		    if(it->get_bit_string().count()>1)
		    {
			iterator final=data.begin()+((it-data.begin())/flavours+1)*flavours;
			do
			{
			    --final;
			}
			while(!(final->is_marked()));
			if(it==final)
			{
			    it=data.begin()+((it-data.begin())/flavours)*flavours;
			    while(!(it->is_marked()))
			    {
				++it;
			    }
			    return true;
			}
			else
			{
			    do
			    {
				++it;
			    }
			    while(!(it->is_marked()));
			    return false;
			}
		    }
		    return true;
		}
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"next current iterator requested for uninitialised data"<<endlog;
		return false;
	    }

	    /* Function finding the next series of marked currents w.r.t. the
	     * argument. It starts performing next_marked_current from the back
	     * of the vector, contuining to the front as long as the iterators
	     * were the last marked currents in their channel block: */

	    static bool next_marked_currents(std::vector<iterator>& iters)
	    {
		iterator check=first_marked_current(iters[0]->get_bit_string());
		size_type n=iters.size()-1;
		while(next_marked_current(iters[n]) and n != 0)
		{
		    --n;
		}
		return (n != 0 or iters[0] != check);
	    }

	    /* Function unmarking all currents: */

	    static void unmark()
	    {
		for(size_type i=0;i<data.size();++i)
		{
		    data[i].unmark();
		}
	    }

	private:

	    /* Vector of currents of all possible flavours in all possible
	     * momentum channels: */

	    static std::vector<current_type>data;
	    
	    /* Number of flavours in the model: */
	    
	    static size_type flavours;

	    /* Initialisation tag: */

	    static bool initialised;
    };

    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t current_tree<model_t,N_in,N_out>::N_bits;
    template<class model_t,std::size_t N_in,std::size_t N_out>const std::size_t current_tree<model_t,N_in,N_out>::static_size;
    template<class model_t,std::size_t N_in,std::size_t N_out>std::size_t current_tree<model_t,N_in,N_out>::N_final=0;
    template<class model_t,std::size_t N_in,std::size_t N_out>std::vector<typename current_tree<model_t,N_in,N_out>::current_type >current_tree<model_t,N_in,N_out>::data;
    template<class model_t,std::size_t N_in,std::size_t N_out>typename current_tree<model_t,N_in,N_out>::size_type current_tree<model_t,N_in,N_out>::flavours=0;
    template<class model_t,std::size_t N_in,std::size_t N_out>bool current_tree<model_t,N_in,N_out>::initialised=false;
}

#endif /*CAMGEN_CURRENT_TREE_H_*/


