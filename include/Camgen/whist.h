//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file whist.h
    \brief Weight histogrammer facility for MC generators.
 */

#ifndef CAMGEN_WHIST_H_
#define CAMGEN_WHIST_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Weight histogrammer facility for Monte Carlo generators. It is used to    *
 * replace the true maximum weight by a reduced wmax which yields histograms *
 * that are identical up to a given accuracy.                                *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <vector>
#include <Camgen/debug.h>
#include <Camgen/logstream.h>
#include <Camgen/plt_config.h>
#include <Camgen/plt_script.h>

namespace Camgen
{
    /// Weight histogramming class for maximal epsilon-weight determination.

    template<class value_t>class weight_histogrammer
    {
	public:

	    /* Type definition: */

	    typedef value_t value_type;
	    typedef typename std::vector<value_t>::size_type size_type;

	    /* Public constructors: */
	    /*----------------------*/

	    /// Constructor with maximum weight and number of bins.

	    weight_histogrammer(const value_type& maxw_,size_type nbins):freqs(std::max((size_type)2,nbins),(value_type)0),maxw(maxw_),wsum(maxw_)
	    {
		freqs.back()=(value_type)1;
	    }

	    /// Constructor with number of bins argument (max weight is taken 1).

	    weight_histogrammer(size_type nbins):freqs(std::max((size_type)2,nbins),(value_type)0),maxw(1),wsum(0){}
	    

	    /* Public modifiers: */
	    /*-------------------*/

	    /// Inserts a weight. If the argument is larger than the current
	    /// histogram allows, rebins the histogram such that the number of bins is
	    /// conserved.
	    
	    void insert(const value_type& w)
	    {
		if(w<(value_type)0 or freqs.size()<2)
		{
		    return;
		}
		value_type bw(bin_width());
		if(w<maxw)
		{
		    ++(freqs[size_type(w/bw)]);
		}
		else
		{
		    rebin(w);
		    maxw=w;
		    ++(freqs.back());
		}
		wsum+=w;
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Number of bins.

	    size_type bins() const
	    {
		return freqs.size();
	    }

	    /// Bin width.

	    value_type bin_width() const
	    {
		if(freqs.size()<2)
		{
		    return 0;
		}
		return maxw/(bins()-1);
	    }

	    /// Minimum value histogrammmed.

	    value_type lower_bound() const
	    {
		return 0;
	    }

	    /// Maximal value histogrammed.

	    value_type upper_bound() const
	    {
		return (maxw+bin_width());
	    }

	    /// Sum of all bin contents.

	    value_type events() const
	    {
		value_type result(0);
		for(size_type i=0;i<freqs.size();++i)
		{
		    result+=freqs[i];
		}
		return result;
	    }

	    /// Sum of all bin contents times bin width.

	    value_type integral() const
	    {
		value_type result(0);
		for(size_type i=0;i<freqs.size();++i)
		{
		    result+=freqs[i]*i;
		}
		return result*bin_width();
	    }

	    /// Maximal frequency.

	    value_type max_freq() const
	    {
		value_type result(0);
		for(size_type i=0;i<freqs.size();++i)
		{
		    result=std::max(result,freqs[i]);
		}
		return result;
	    }

	    /// Computes the weight for which all higher-weight events make up
	    /// epsilon times the total cross section.

	    value_type max_weight(const value_type& epsilon) const
	    {
		if(!(epsilon>(value_type)0))
		{
		    return maxw;
		}
		if(!(epsilon<(value_type)1))
		{
		    return 0;
		}
		typename std::vector<value_type>::const_reverse_iterator it=freqs.rbegin();
		size_type i=freqs.size()-1;
		value_type bw(bin_width()),bound(epsilon*wsum/bw),Splus((*it)*i);
		while(Splus<bound and it!=freqs.rend())
		{
		    ++it;
		    --i;
		    Splus+=((*it)*i);
		}
		return (i-1+(Splus-bound)/(i*(*it)))*bw;
	    }

	    /// Plot script.

	    plot_script* plot(const std::string& filename,const char* term=NULL) const
	    {
		plot_script* result=new plot_script(filename,(double)0,static_cast<double>(maxw),(double)0.1,(double)1.5*static_cast<double>(max_freq()),term);
		result->histo=true;
		data_wrapper* data;
		value_type bin,val;
		if(plot_config::gnuplot_path==NULL)
		{
		    data=new data_wrapper(filename+".dat",&bin,&val);
		}
		else
		{
		    data=new data_wrapper(&bin,&val);
		}
		value_type bw(bin_width());
		for(size_type i=0;i<freqs.size();++i)
		{
		    bin=i*bw;
		    val=freqs[i]*bw;
		    data->fill();
		}
		data->write();
		data_stream* dstr=new data_stream(data);
		dstr->title="weights";
		dstr->style="histeps";
		result->add_plot(dstr);
		result->ylog=true;
		return result;
	    }

	protected:

	    /* Protected constructors: */
	    /*-------------------------*/

	    /// Default constructor.

	    weight_histogrammer():maxw(1),wsum(0){}

	private:

	    /* Rebinning utility class. */

	    class bin_iterator
	    {
		public:

		    typedef typename std::vector<value_type>::iterator iterator;

		    bin_iterator(iterator pos_):pos(pos_),offset(0){}

		    value_type integrate(value_type range,iterator end)
		    {
			if(pos==end)
			{
			    return 0;
			}
			if(range+offset<1)
			{
			    offset+=range;
			    return range*(*pos);
			}
			value_type r=range+offset-1;
			value_type result=(1-offset)*(*pos);
			size_type n=static_cast<size_type>(std::floor(r));
			for(size_type i=0;i<n;++i)
			{
			    if((++pos)==end)
			    {
				break;
			    }
			    result+=(*pos);
			}
			if(pos!=end)
			{
			    offset=r-n;
			    ++pos;
			    result+=(offset*(*pos));
			}
			else
			{
			    offset=0;
			}
			return result;
		    }
		
		    iterator pos;
		    value_type offset;
	    };

	    /* Rebins such that the new maximal weight w is contained in the
	     * last bin and the number of bins is not changed. */

	    void rebin(const value_type& w)
	    {
		if(w<maxw)
		{
		    return;
		}
		value_type oldwidth=bin_width();
		value_type newwidth=w/(freqs.size()-1);
		std::vector<value_type>newfreqs(freqs.size());
		bin_iterator it(freqs.begin());
		for(size_type i=0;i<newfreqs.size()-1;++i)
		{
		    newfreqs[i]=it.integrate(newwidth/oldwidth,freqs.end());
		}
		value_type lastfill=it.integrate(newwidth/oldwidth,freqs.end());
		if(lastfill>(value_type)0)
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"the last bin is incorrectly filled with "<<lastfill<<endlog;
		}
		freqs=newfreqs;
	    }

	    /* Histogram content: */

	    std::vector<value_type>freqs;

	    /* Maximal weight=final bin position: */

	    value_type maxw;

	    /* Weight sum: */

	    value_type wsum;
    };
}

#endif /*CAMGEN_WHIST_H_*/
