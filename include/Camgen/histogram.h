//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file histogram.h
    \brief Small histogramming class definition.
 */

#ifndef CAMGEN_HISTOGRAM_H_
#define CAMGEN_HISTOGRAM_H_

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <Camgen/plt_config.h>
#include <Camgen/plt_script.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Histogramming facility for Camgen's random number generators. The class      *
 * performs a binning of the data in a range offset*sigma around the mean value, *
 * over a given number of bins. Then it writes out the bin fillings and the bin  *
 * weights as given by the generator to a data file, together with a gnuplot     *
 * script to visualise the histogram.                                            *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Histogram-producing class template. The first template argument denotes the type of
    /// the random number, the second argument a comparison class.

    template<class T,class comp_type=std::less<T> >class histogram
    {
	public:

	    /* Useful type definitions: */

	    typedef T value_type;
	    typedef comp_type comparison_type;
	    
	    /// Data point type definition.

	    class data_point
	    {
		public:
		    data_point(const value_type& val_,const value_type wght_=(value_type)1):val(val_),wght(wght_){}
		    value_type& value()
		    {
			return val;
		    }
		    const value_type& value() const
		    {
			return val;
		    }
		    value_type& weight()
		    {
			return wght;
		    }
		    const value_type& weight() const
		    {
			return wght;
		    }
		private:
		    value_type val,wght;
	    };

	    /// Data point comparison function according to value.

	    static bool data_value_comp(const data_point& x1,const data_point& x2)
	    {
		return x1.value()<x2.value();
	    }
	    
	    /// Data point comparison function according to weight.

	    static bool data_weight_comp(const data_point& x1,const data_point& x2)
	    {
		return x1.weight()<x2.weight();
	    }

	    /// Bin class type definition.

	    class bin
	    {
		public:
		    bin(const value_type& position_,std::size_t filling_=0,const value_type& wght_=(value_type)1):position(position_),filling(filling_),wght(wght_){}
		    value_type& pos()
		    {
			return position;
		    }
		    const value_type& pos() const
		    {
			return position;
		    }
		    std::size_t& value()
		    {
			return filling;
		    }
		    const std::size_t& value() const
		    {
			return filling;
		    }
		    value_type& weight()
		    {
			return wght;
		    }
		    const value_type& weight() const
		    {
			return wght;
		    }
		private:
		    value_type position;
		    std::size_t filling;
		    value_type wght;
	    };

	    /// Bin comparison by position.

	    static bool bin_pos_comp(const bin& b1,const bin& b2)
	    {
		return b1.pos()<b2.pos();
	    }
	    
	    /// Bin postion-value comparison.

	    static bool value_pos_comp(const value_type& x,const bin& b)
	    {
		return x<b.pos();
	    }

	    /// Bin comparison by value.

	    static bool bin_value_comp(const bin& b1,const bin& b2)
	    {
		return b1.value()<b2.value();
	    }

	    /// Bin comparison by weight.

	    static bool bin_weight_comp(const bin& b1,const bin& b2)
	    {
		return b1.weight()<b2.weight();
	    }

	    typedef typename std::vector<data_point>::size_type data_size_type;
	    typedef typename std::vector<data_point>::iterator data_iterator;
	    typedef typename std::vector<data_point>::reverse_iterator reverse_data_iterator;
	    typedef typename std::vector<data_point>::const_iterator const_data_iterator;
	    typedef typename std::vector<data_point>::const_reverse_iterator const_reverse_data_iterator;
	    typedef typename std::vector<bin>::size_type bin_size_type;
	    typedef typename std::vector<bin>::iterator bin_iterator;
	    typedef typename std::vector<bin>::reverse_iterator reverse_bin_iterator;
	    typedef typename std::vector<bin>::const_iterator const_bin_iterator;
	    typedef typename std::vector<bin>::const_reverse_iterator const_reverse_bin_iterator;

	    /// Constructor: the first argument is a pointer whose value will be
	    /// written to the histogram each time the store-function is called.
	    /// The second argument is optional and provides a speedup for large
	    /// samples, when taken to be the expected sample size.

	    histogram(const value_type* x,const value_type* w=NULL,std::size_t pts=1000):var(x),weight(w),N_points(0),mu(0),weight_sum(0)
	    {
		data.reserve(pts);
	    }

	    /// Stores the pointer value in the histogram.

	    void store()
	    {
		if(weight==NULL)
		{
		    data.push_back(data_point(*var));
		    ++weight_sum;
		}
		else
		{
		    data.push_back(data_point(*var,*weight));
		    weight_sum+=(*weight);
		}
		mu+=(*var);
		++N_points;
	    }

	    /// Creates a histogram from the data. The first argument determines
	    /// the number of bins to be used, the second and third resp. the
	    /// number of standard deviations from left and right of the mean to
	    /// be plotted.

	    void make(bin_size_type bns=100,const value_type& offset_left=3,const value_type& offset_right=3)
	    {
		N_bins=bns;
		bins.clear();
		bins.reserve(N_bins);
		mu/=(N_points);
		sig=0;
		for(data_size_type i=0;i<N_points;++i)
		{
		    sig+=((data[i].value()-mu)*(data[i].value()-mu));
		}
		sig=std::sqrt(sig/N_points);
		
		/* Computing minimal and maximal bin and bin size: */
		
		x_min=mu-std::abs(offset_left)*sig;
		x_max=mu+std::abs(offset_right)*sig;
		bin_size=(x_max-x_min)/N_bins;

		/* Constructing the bins: */

		if(weight==NULL)
		{
		    for(bin_size_type i=0;i<N_bins;++i)
		    {
			bins.push_back(bin(x_min+i*bin_size));
		    }
		}
		else
		{
		    for(bin_size_type i=0;i<N_bins;++i)
		    {
			bins.push_back(bin(x_min+i*bin_size,0,0));
		    }
		}

		/* Filling bins: */

		if(weight==NULL)
		{
		    for(data_size_type i=0;i<N_points;++i)
		    {
			bin b(data[i].value());
			typename std::vector<bin>::iterator it=std::lower_bound(bins.begin(),bins.end(),b,bin_pos_comp);
			if(it != bins.begin() and it != bins.end())
			{
			    ++it->value();
			}
		    }
		}
		else
		{
		    for(data_size_type i=0;i<N_points;++i)
		    {
			bin b(data[i].value(),0,0);
			typename std::vector<bin>::iterator it=std::lower_bound(bins.begin(),bins.end(),b,bin_pos_comp);
			if(it != bins.begin() and it != bins.end())
			{
			    ++it->value();
			    it->weight()+=data[i].weight();
			}
		    }
		    value_type norm=N_points/weight_sum;
		    for(bin_size_type i=0;i<N_bins;++i)
		    {
			bins[i].weight()*=norm;
		    }
		}
		max_freq=std::max_element(bins.begin(),bins.end(),bin_value_comp)->value();
	    }

	    /// Outputs the bin positions, values and weights to the output
	    /// stream argument.

	    std::ostream& write(std::ostream& os) const
	    {
		if(weight==NULL)
		{
		    for(bin_size_type i=0;i<N_bins;++i)
		    {
			os<<bins[i].pos()<<"\t"<<bins[i].value()<<std::endl;
		    }
		}
		else
		{
		    for(bin_size_type i=0;i<N_bins;++i)
		    {
			os<<bins[i].pos()<<"\t"<<bins[i].value()<<"\t"<<bins[i].weight()<<std::endl;
		    }
		}
		return os;
	    }

	    /// Writes the gnuplot script for the histogram plot.

	    plot_script* plot(const std::string& filename,const char* term=NULL) const
	    {
		plot_script* result=new plot_script(filename,(double)x_min,(double)x_max,(double)0,(double)10*(double)max_freq,term);
		result->histo=true;
		data_wrapper* data;
		value_type bin,val,wgt;
		if(weight==NULL)
		{
		    if(plot_config::gnuplot_path==NULL)
		    {
			data=new data_wrapper(filename+".dat",&bin,&val);
		    }
		    else
		    {
			data=new data_wrapper(&bin,&val);
		    }
		    for(bin_size_type i=0;i<N_bins;++i)
		    {
			bin=bins[i].pos();
			val=bins[i].value();
			data->fill();
		    }
		    data->write();
		    data_stream* dstr=new data_stream(data);
		    dstr->title="generated";
		    dstr->style="histeps";
		    result->add_plot(dstr);
		}
		else
		{
		    if(plot_config::gnuplot_path==NULL)
		    {
			data=new data_wrapper(filename+".dat",&bin,&val,&wgt);
		    }
		    else
		    {
			data=new data_wrapper(&bin,&val,&wgt);
		    }
		    for(bin_size_type i=0;i<N_bins;++i)
		    {
			bin=bins[i].pos();
			val=bins[i].value();
			wgt=bins[i].weight();
			data->fill();
		    }
		    data->write();
		    data_stream* dstr=new data_stream(data,"2");
		    dstr->title="generated";
		    dstr->style="histeps";
		    result->add_plot(dstr);
		    data_stream* wstr=new data_stream(data,"3");
		    wstr->title="weights";
		    wstr->style="histeps";
		    result->add_plot(wstr);
		}
		return result;
	    }

	    /// Returns the mean.

	    value_type mean() const
	    {
		return mu;
	    }

	    /// Returns the variance.

	    value_type sigma() const
	    {
		return sig;
	    }

	    /// Returns the number of data points.

	    std::size_t N_events() const
	    {
		return N_points;
	    }

	    /// Const iterator to the first data point.

	    const_data_iterator data_begin() const
	    {
		return data.begin();
	    }

	    /// Const iterator past the last data point.

	    const_data_iterator data_end() const
	    {
		return data.end();
	    }

	    /// Const reverse iterator to the last data point.

	    const_reverse_data_iterator data_rbegin() const
	    {
		return data.rbegin();
	    }

	    /// Const reverse iterator before the first data point.

	    const_reverse_data_iterator data_rend() const
	    {
		return data.rend();
	    }

	    /// Const iterator to the first bin.

	    const_bin_iterator bins_begin() const
	    {
		return bins.begin();
	    }

	    /// Const iterator past the last bin.

	    const_bin_iterator bins_end() const
	    {
		return bins.end();
	    }

	    /// Const reverse iterator to the last bin.

	    const_reverse_bin_iterator bins_rbegin() const
	    {
		return bins.rbegin();
	    }

	    /// Const reverse iterator before the first bin.

	    const_reverse_bin_iterator bins_rend() const
	    {
		return bins.rend();
	    }

	    /// Returns the bin size.

	    value_type binsize() const
	    {
		return bin_size;
	    }

	    /// Counts events in the histogram.

	    data_size_type count_events() const
	    {
		data_size_type sum(0);
		for(const_bin_iterator it=bins.begin();it!=bins.end();++it)
		{
		    sum+=it->value();
		}
		return sum;
	    }

	private:

	    /* Variable for which the distribution is computed: */

	    const value_type* var;

	    /* Weight of the stochastic variable value: */

	    const value_type* weight;

	    /* Number of Monte Carlo points: */

	    data_size_type N_points;

	    /* Number of bins: */

	    bin_size_type N_bins;
	    
	    /* Upper and lower x-values: */
	    
	    value_type x_min,x_max;

	    /* Mean and variance and weight sum of the point set: */

	    value_type mu,sig,weight_sum;
	    
	    /* Width of the bins: */
	    
	    value_type bin_size;
	    
	    /* Stored data points: */

	    std::vector<data_point>data;

	    /* Bin x-values: */

	    std::vector<bin>bins;
	    
	    /* Maximum of the histogram: */
	    
	    value_type max_freq;

	    /* Comparison object: */

	    comp_type less_than;
    };
}

#endif /*CAMGEN_HISTOGRAM_H_*/


