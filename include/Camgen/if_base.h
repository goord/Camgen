//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file if_base.h
  \brief Abstract base class for output interfaces.
  */

#ifndef CAMGEN_IF_BASE_H_
#define CAMGEN_IF_BASE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Abstract base class for event generator interface classes. The base class   *
 * contains the abstract fill()/write()/reset() methods, plus some common data *
 * members such as the number of streamed events.                              *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <cstdlib>
#include <vector>
#include <Camgen/ps_gen_base.h>

namespace Camgen
{
    /// Abstract base class for interfaces.

    template<class model_t>class interface_base
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef ps_generator_base<model_t> generator_type;
	    typedef typename generator_type::momentum_type momentum_type;
	    typedef typename generator_type::value_type value_type;
	    typedef typename generator_type::size_type size_type;

	    struct stats_entry
	    {
		size_type input_evts;
		size_type output_evts;
		MC_integral<value_type>xsec;

		stats_entry():input_evts(0),output_evts(0){}
	    };

	    /* Constructors/destructors: */
	    /*---------------------------*/

	    /// Constructor.

	    interface_base(generator_type* gen_):gen(gen_),input_evts(0),output_evts(0),zero_weight_flag(false){}

	    /// Destructor.

	    virtual ~interface_base(){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Abstract event streaming method.

	    virtual bool fill()
	    {
		proc_id=gen->process_id();
		++input_evts;
		stats_entry& stats=proc_stats[proc_id];
		++(stats.input_evts);
		stats.xsec=gen->sub_xsec();
		w=gen->w();
		if(w<(value_type)0)
		{
		    return false;
		}
		if(w==(value_type)0 and !zero_weight_flag)
		{
		    return false;
		}
		if(fill_event())
		{
		    ++output_evts;
		    ++stats.output_evts;
		    return true;
		}
		return false;
	    }

	    /// Abstract output file writing method.

	    virtual bool write()=0;

	    /// Sets the zero-weight flag.

	    virtual void write_zero_weight_events()
	    {
		zero_weight_flag=false;
	    }

	    /// Unsets the zero-weight flag.

	    virtual void skip_zero_weight_events()
	    {
		zero_weight_flag=false;
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Returns the number of subprocesses.

	    size_type processes() const
	    {
		return (gen==NULL)?0:(gen->processes());
	    }

	    /// Returns the event size.

	    virtual size_type event_size() const=0;

	    /// Returns the number of events streamed to the interface.

	    size_type input_events() const
	    {
		return input_evts;
	    }

	    /// Returns the number of events streamed to the interface for process i.

	    size_type input_events(int i) const
	    {
		typename std::map<int,stats_entry>::const_iterator it=proc_stats.find(i);
		return (it==proc_stats.end())?0:(it->second.input_evts);
	    }

	    /// Returns the number of events written to disk.

	    size_type output_events() const
	    {
		return output_evts;
	    }

	    /// Returns the number of events written to disk for process i.

	    size_type output_events(int i) const
	    {
		typename std::map<int,stats_entry>::const_iterator it=proc_stats.find(i);
		return (it==proc_stats.end())?0:(it->second.output_evts);
	    }

	    /// Returns the total size of the output file.
	    
	    size_type file_size() const
	    {
		return output_evts*event_size();
	    }

	    /// Returns whether we are writing zero-weight events.

	    bool writes_zero_weights() const
	    {
		return zero_weight_flag;
	    }

	    /// Writes the statistics file, containing a list with number of
	    /// generated events, cross section and error per subprocess.

	    std::ostream& print_statistics(std::ostream& os) const
	    {
		for(typename std::map<int,stats_entry>::const_iterator it=proc_stats.begin();it!=proc_stats.end();++it)
		{
		    os<<std::setw(20)<<std::left<<it->first;
		    os<<std::setw(20)<<std::left<<it->second.input_evts;
		    os<<std::setw(20)<<std::left<<it->second.xsec.value;
		    os<<std::setw(20)<<std::left<<it->second.xsec.error<<std::endl;
		}
		return os;
	    }

	protected:

	    /* Generator instance: */

	    generator_type* gen;

	    /* Event weight: */

	    value_type w;

	    /* Subprocess id: */

	    int proc_id;

	    /* Abstract method filling the event variables: */

	    virtual bool fill_event()=0;

	private:

	    /* Number of events streamed to the interface. */

	    size_type input_evts;

	    /* Number of events written to disk. */

	    size_type output_evts;

	    /* Flag denoting whether events with zero weight are written to tape. */

	    bool zero_weight_flag;

	    /* Map with process-id keys and stats entries as values. */

	    std::map<int,stats_entry> proc_stats;
    };
}

#endif /*CAMGEN_IF_BASE_H_*/

