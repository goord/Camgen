//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file qcd_cols.h
    \brief QCD colour generators class interface and implementation.
 */

#ifndef CAMGEN_QCD_COLS_H_
#define CAMGEN_QCD_COLS_H_

#include <Camgen/adjoint.h>
#include <Camgen/col_flow.h>
#include <Camgen/col_gen.h>
#include <Camgen/uni_sphere.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Colour generators for QCD-like models, which have one unbroken gauge group    *
 * and particles in either the fundamental, adjoint or singlet representation.   *
 * Such models are required to contain a static constant parameter 'N_c',        *
 * denoting the dimension of the fundamental representation. There are 2         *
 * classes, one for the gluons in either the adjoint or colour flow              *
 * representation. Each one is specialised for continuous and discrete colours.  *
 * The colour-flow, discrete colours generator returns only configurations which *
 * preserve colour and is therefore not uniform.                                 *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// Uniform colour generator for QCD-like theories with gluons in the
    /// adjoint representation. The first template parameter denotes the model
    /// type, the second and third integers respectively the initial and final
    /// state multiplicities, the fourth integer the dimension of the
    /// fundamental representation, the fifth parameter the random number
    /// generator and the last boolean the continuous colours tag.

    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t,bool cont_colours>class adjoint_QCD;

    /// Adjoint-QCD colour generator specialisation for discrete colours.

    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>class adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,false>: public colour_generator<value_t,N_in,N_out,false>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/
	    
	    typedef value_t value_type;
	    typedef colour_generator<value_t,N_in,N_out,false> base_type;
	    typedef typename base_type::size_type size_type;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    
	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Factory method: */

	    template<class model_t>static adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,false>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,false>* result=new adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,false>(base_type::template make_colour_vector<model_t>(it),base_type::template make_colour_range_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }
	    
	    /* Utility function: */
	    
	    static vector<std::vector<size_type>,N_tot> make_rank_vector(vector<int,N_tot>types)
	    {
		vector<std::vector<size_type>,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    if(types[i]==1 or types[i]==-1)
		    {
			result[i]=std::vector<size_type>(1,N_c);
		    }
		    if(types[i]==2 or types[i]==-2)
		    {
			result[i]=std::vector<size_type>(1,N_c*N_c-1);
		    }
		}
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/
	    
	    /* Non allocating constructor: */

	    adjoint_QCD(const vector<std::vector<size_type>*,N_tot>& cols_,const vector<std::vector<size_type>,N_tot>& ranges_):base_type(cols_,ranges_),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<this->colour_rank(i);++j)
		    {
			const_weight*=(this->colour_range(i,j));
		    }
		}
	    }

	    /* Non allocating constructor. The second argument is an integer
	     * vector denoting the representation types of the external
	     * particles: 0 for singlets, 1 or -1 for fundamental representations
	     * and 2 or -2 for adjoint representations. */

	    adjoint_QCD(const vector<std::vector<size_type>*,N_tot>& cols_,const vector<int,N_tot>& types):base_type(cols_,make_rank_vector(types)),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<this->colour_rank(i);++j)
		    {
			const_weight*=(this->colour_range(i,j));
		    }
		}
	    }

	    /* Allocating constructor. The argument is an integer
	     * vector denoting the representation types of the external
	     * particles: 0 for singlets, 1 or -1 for fundamental representations
	     * and 2 or -2 for adjoint representations. */

	    adjoint_QCD(const vector<int,N_tot>& types):base_type(make_rank_vector(types)),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<this->colour_rank(i);++j)
		    {
			const_weight*=(this->colour_range(i,j));
		    }
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generation function: */

	    bool generate()
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    for(size_type j=0;j<this->colour_rank(i);++j)
		    {
			this->colour(i,j)=rn_stream::throw_dice(this->colour_range(i,j));
		    }
		}
		this->weight()=(value_type)const_weight;
		return true;
	    }

	    /* Weight evaluation function: */

	    bool evaluate_weight()
	    {
		this->weight()=(value_type)const_weight;
		return true;
	    }

	    /* Public const methods: */
	    /*----------------------*/

	    /* Clone method: */

	    adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,false>* clone() const
	    {
		return new adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,false>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "qcd";
	    }

	private:

	    size_type const_weight;
    };
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const std::size_t adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,false>::N_tot;

    /* Adjoint-QCD colour generator specialisation for continuous colours: */

    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>class adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>: public colour_generator<value_t,N_in,N_out,true>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/

	    typedef value_t value_type;
	    typedef colour_generator<value_t,N_in,N_out,true> base_type;
	    typedef std::complex<value_t> c_value_type;
	    typedef typename base_type::object_type object_type;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename base_type::size_type size_type;
	    
	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Factory method: */
	    
	    template<class model_t>static adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>* result=new adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>(base_type::template make_colour_tensor_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }

	    /* Utility function: */
	    
	    static vector<std::vector<size_type>,N_tot> make_rank_vector(vector<int,N_tot>types)
	    {
		vector<std::vector<size_type>,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    if(types[i]==1 or types[i]==-1)
		    {
			result[i]=std::vector<size_type>(1,N_c);
		    }
		    if(types[i]==2 or types[i]==-2)
		    {
			result[i]=std::vector<size_type>(1,N_c*N_c-1);
		    }
		}
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Non-allocating standard constructor: */
	    
	    adjoint_QCD(const vector<object_type*,N_tot>& cols_):base_type(cols_){}
	    
	    /* Non-allocating constructor. The second argument is an integer
	     * vector denoting the representation types of the external
	     * particles: 0 for singlets, 1 or -1 for fundamental representations
	     * and 2 or -2 for adjoint representations: */

	    adjoint_QCD(const vector<object_type*,N_tot>& cols_,const vector<int,N_tot>& types):base_type(cols_,make_rank_vector(types)){}

	    /* Allocating constructor. The argument is an integer
	     * vector denoting the representation types of the external
	     * particles: 0 for singlets, 1 or -1 for fundamental representations
	     * and 2 or -2 for adjoint representations: */

	    adjoint_QCD(const vector<int,N_tot>& types):base_type(make_rank_vector(types)){}

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generation operator implementation. Generates fundamental
	     * representation tensors by uniformly distributed points on the
	     * complex N-dimensional sphere and adjoint representations by
	     * uniformly distributed points on the real (N^2-1)-dimensional
	     * sphere: */

	    bool generate()
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    if(this->colour_range(i,0)==N_c)
		    {
			q_gen.generate();
			for(size_type j=0;j<N_c;++j)
			{
			    this->colour_entry(i,j)=q_norm*std::complex<value_type>(qvec[j<<1],qvec[(j<<1)+1]);
			}
		    }
		    if(this->colour_range(i,0)==N_c*N_c-1)
		    {
			g_gen.generate();
			for(size_type j=0;j<N_c*N_c-1;++j)
			{
			    this->colour_entry(i,j)=std::complex<value_type>(g_norm*gvec[j],0);
			}
		    }
		}
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Weight evaluation method. Assigns unity to the weight. */

	    bool evaluate_weight()
	    {
		this->weight()=(value_type)1;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method: */

	    adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>* clone() const
	    {
		return new adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    std::string type() const
	    {
		return "qcd";
	    }

	private:

	    /* Private static data: */
	    /*----------------------*/

	    /* Utility vectors: */

	    static vector<value_t,2*N_c>qvec;
	    static vector<value_t,N_c*N_c-1>gvec;

	    /* Helper generators: */

	    static uniform_sphere<value_t,2*N_c-1,rng_t>q_gen;
	    static uniform_sphere<value_t,N_c*N_c-2,rng_t>g_gen;

	    /* Normalisation factors: */

	    static const value_t q_norm;
	    static const value_t g_norm;
    };
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const std::size_t adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::N_tot;
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>vector<value_t,2*N_c> adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::qvec;
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>vector<value_t,N_c*N_c-1> adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::gvec;
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>uniform_sphere<value_t,2*N_c-1,rng_t> adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::q_gen(&adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::qvec);
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>uniform_sphere<value_t,N_c*N_c-2,rng_t> adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::g_gen(&adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::gvec);
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const value_t adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::q_norm=std::sqrt((value_t)N_c);
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const value_t adjoint_QCD<value_t,N_in,N_out,N_c,rng_t,true>::g_norm=std::sqrt(value_t(N_c*N_c-1));
    
    /// Colour generator for QCD-like threories with gluons in the colour-flow
    /// representation. The first template parameter denotes the model type, the
    /// second one the model type, the second and third integers resp. the initial
    /// and final state multiplicities, the fourth parameter the ranom number
    /// generator and the final boolean the continuous_colours flag.

    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t,bool continuous_colours>class colour_flow_QCD;

    /* Colour-flow QCD-like colour generator specialisation for discrete
    * colours. */
    
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>class colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,false>: public colour_generator<value_t,N_in,N_out,false>
    {
	public:

	    /* Useful type definitions: */
	    /*--------------------------*/
	    
	    typedef value_t value_type;
	    typedef colour_generator<value_t,N_in,N_out,false> base_type;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename base_type::size_type size_type;

	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static functions: */
	    /*--------------------------*/
	    
	    /* Utility functions: */
	    
	    static vector<size_type,N_tot> make_rank_vector(vector<int,N_tot>types)
	    {
		vector<size_type,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    if(types[i]==1 or types[i]==-1)
		    {
			result[i]=1;
		    }
		    else if(types[i]==2 or types[i]==-2)
		    {
			result[i]=2;
		    }
		    else
		    {
			result[i]=0;
		    }
		}
		return result;
	    }

	    template<class model_t>static vector<std::vector<int>,N_tot> make_type_vector(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		vector<std::vector<int>,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    result[i].resize(it->get_phase_space(i)->colour_rank());
		    for(size_type j=0;j<result[i].size();++j)
		    {
			result[i][j]=it->get_phase_space(i)->particle_type->colour_type(j);
		    }
		}
		return result;
	    }

	    /* Factory method: */

	    template<class model_t>static colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,false>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,false>* result=new colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,false>(base_type::template make_colour_vector<model_t>(it),make_type_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/

	    /* Standalone-mode constructor. Argument array entry 0 is interpreted as a
	     * singlet, 1 or -1 as fundamental representation particles, and 2 for adjoint
	     * representations. */

	    colour_flow_QCD(vector<int,N_tot>types):base_type(N_c,make_rank_vector(types))
	    {
		size_type wf_denom(1);
		for(size_type i=0;i<N_in;++i)
		{
		    switch(types[i])
		    {
			case 1:
			    colours.push_back(std::pair<size_type,int>(i,0));
			    break;
			case -1:
			    anti_colours.push_back(std::pair<size_type,int>(i,0));
			    break;
			case 2:
			    colours.push_back(std::pair<size_type,int>(i,0));
			    anti_colours.push_back(std::pair<size_type,int>(i,1));
			    wf_denom*=2;
			    break;
			case -2:
			    colours.push_back(std::pair<size_type,int>(i,0));
			    anti_colours.push_back(std::pair<size_type,int>(i,1));
			    wf_denom*=2;
			    break;
		    }
		}
		for(size_type i=N_in;i<N_tot;++i)
		{
		    switch(types[i])
		    {
			case 1:
			    anti_colours.push_back(std::pair<size_type,int>(i,0));
			    break;
			case -1:
			    colours.push_back(std::pair<size_type,int>(i,0));
			    break;
			case 2:
			    colours.push_back(std::pair<size_type,int>(i,1));
			    anti_colours.push_back(std::pair<size_type,int>(i,0));
			    wf_denom*=2;
			    break;
			case -2:
			    colours.push_back(std::pair<size_type,int>(i,1));
			    anti_colours.push_back(std::pair<size_type,int>(i,0));
			    wf_denom*=2;
			    break;
		    }
		}
		if(colours.size()!=anti_colours.size())
		{
		    log(log_level::warning)<<CAMGEN_STREAMLOC<<"process allows no colour-conserving flows: all returned weight will be 0"<<endlog;
		    weight_factor=(value_type)0;
		}
		else
		{
		    size_type wf_num(1);
		    for(size_type i=0;i<colours.size();++i)
		    {
			wf_num*=(N_c*(i+1));
		    }
		    weight_factor=(value_type)wf_num/(value_type)wf_denom;
		}
	    }

	    /* Non-allocating constructor. */

	    colour_flow_QCD(const vector<std::vector<size_type>*,N_tot>& cols_,const vector<std::vector<int>,N_tot>& types):base_type(cols_,N_c)
	    {
		size_type wf_denom(1);
		for(size_type i=0;i<N_in;++i)
		{
		    if(this->colour_rank(i)==1)
		    {
			if(types[i][0]==1)
			{
			    colours.push_back(std::pair<size_type,int>(i,0));
			}
			if(types[i][0]==-1)
			{
			    anti_colours.push_back(std::pair<size_type,int>(i,0));
			}
		    }
		    else if(this->colour_rank(i)==2)
		    {
			colours.push_back(std::pair<size_type,int>(i,0));
			anti_colours.push_back(std::pair<size_type,int>(i,1));
			wf_denom*=2;
		    }
		}
		for(size_type i=N_in;i<N_tot;++i)
		{
		    if(this->colour_rank(i)==1)
		    {
			if(types[i][0]==1)
			{
			    anti_colours.push_back(std::pair<size_type,int>(i,0));
			}
			if(types[i][0]==-1)
			{
			    colours.push_back(std::pair<size_type,int>(i,0));
			}
		    }
		    else if(this->colour_rank(i)==2)
		    {
			anti_colours.push_back(std::pair<size_type,int>(i,0));
			colours.push_back(std::pair<size_type,int>(i,1));
			wf_denom*=2;
		    }
		}
		if(colours.size()!=anti_colours.size())
		{
		    log(log_level::warning)<<"process allows no colour-conserving flows: all returned weight will be 0"<<endlog;
		    weight_factor=(value_type)0;
		}
		else
		{
		    size_type wf_num(1);
		    for(size_type i=0;i<colours.size();++i)
		    {
			wf_num*=(N_c*(i+1));
		    }
		    weight_factor=(value_type)wf_num/(value_type)wf_denom;
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generation operator implementation. This method will only
	     * generate colour-preserving configurations that give (mostly)
	     * nonzero matrix elements in QCD-like models. The weight can vary, 
	     * as this distribution is nonuniform: */

	    bool generate()
	    {
		if(weight_factor==(value_type)0)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		/* Generate the colours randomly: */

		for(size_type i=0;i<N_c;++i)
		{
		    colarray[i]=0;
		}
		size_type c;
		for(size_type i=0;i<colours.size();++i)
		{
		    c=rn_stream::throw_dice(N_c);
		    ++colarray[c];
		    this->colour(colours[i].first,colours[i].second)=c;
		}

		/* Shuffle the anticolour entries: */

		shuffle_partons();
		
		/* Assign the anticolours: */
		
		for(size_type i=0;i<colours.size();++i)
		{
		    this->colour(shuffled_anti_colours[i].first,shuffled_anti_colours[i].second)=this->colour(colours[i].first,colours[i].second);
		}

		/* Compute the corresponding weight: */

		size_type denom=1;
		for(size_type i=0;i<N_c;++i)
		{
		    denom*=factorial(colarray[i]);
		}
		this->weight()=weight_factor/(value_type)denom;
		return true;
	    }

	    /* Event weight evaluation implementation. */

	    bool evaluate_weight()
	    {
		if(weight_factor==(value_type)0)
		{
		    this->weight()=(value_type)0;
		    return false;
		}
		for(size_type n=0;n<N_c;++n)
		{
		    colarray[n]=0;
		}
		for(size_type i=0;i<colours.size();++i)
		{
		    ++colarray[this->colour(colours[i].first,colours[i].second)];
		}
		size_type denom=1;
		for(size_type i=0;i<N_c;++i)
		{
		    denom*=factorial(colarray[i]);
		}
		this->weight()=weight_factor/(value_type)denom;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method: */

	    colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,false>* clone() const
	    {
		return new colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,false>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Les Houches event file output. This will convert the
	     * configuration to an infinite-colour output, starting from colour
	     * number 501: */

	    void LH_output(vector<int,N_tot>& cols,vector<int,N_tot>& anticols) const
	    {
		int n=501;
		cols.assign(0);
		anticols.assign(0);
		for(size_type i=0;i<colours.size();++i)
		{
		    if(colours[i].first<2)
		    {
			cols[colours[i].first]=n;
		    }
		    else
		    {
			anticols[colours[i].first]=n;
		    }
		    if(shuffled_anti_colours[i].first<2)
		    {
			anticols[shuffled_anti_colours[i].first]=n;
		    }
		    else
		    {
			cols[shuffled_anti_colours[i].first]=n;
		    }
		    ++n;
		}
	    }

	    /* Polymorphic type identifier: */

	    std::string type() const
	    {
		return "qcd";
	    }

	private:

	    /* Private data: */
	    /*---------------*/
	    
	    /* Vector of particle numbers that carry colour: */
	    
	    std::vector< std::pair<size_type,int> >colours;

	    /* Vector of particle numbers that carry anticolour: */

	    std::vector< std::pair<size_type,int> >anti_colours;
	    
	    /* Permuted vector of anticoloured particles: */
	    
	    std::vector< std::pair<size_type,int> >shuffled_anti_colours;

	    /* Counter array for colours: */

	    size_type colarray[N_c];

	    /* weight numerator: */

	    value_type weight_factor;

	    /* Private modifiers: */
	    /*--------------------*/

	    /* Function randomly permuting the anticolour entries: */

	    void shuffle_partons()
	    {
		if(anti_colours.size()>0)
		{
		    shuffled_anti_colours=anti_colours;
		    for(size_type i=shuffled_anti_colours.size()-1;i!=0;--i)
		    {
			std::swap(shuffled_anti_colours[i],shuffled_anti_colours[rn_stream::throw_dice(i+1)]);
		    }
		}
	    }
    };
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const std::size_t colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,false>::N_tot;
    
    /* Colour-flow QCD-like colour generator specialisation for continuous
    colours. */

    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>class colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>: public colour_generator<value_t,N_in,N_out,true>
    {
	public:

	    /* Type definitions: */
	    /*-------------------*/
	    
	    typedef value_t value_type;
	    typedef std::complex<value_t> c_value_type;
	    typedef colour_generator<value_t,N_in,N_out,true> base_type;
	    typedef typename base_type::object_type object_type;
	    typedef random_number_stream<value_type,rng_t> rn_stream;
	    typedef typename base_type::size_type size_type;

	    /* Public static data: */
	    /*---------------------*/

	    static const std::size_t N_tot=N_in+N_out;

	    /* Public static methods: */
	    /*------------------------*/

	    /* Factory method: */

	    template<class model_t>static colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>* create_instance(typename CM_algorithm<model_t,N_in,N_out>::tree_iterator it)
	    {
		colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>* result=new colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>(base_type::template make_colour_tensor_vector<model_t>(it));
		result->prefactor=base_type::template make_prefactor<model_t>(it);
		return result;
	    }

	    /* Utility function: */
	    
	    static vector<size_type,N_tot> make_rank_vector(vector<int,N_tot>types)
	    {
		vector<size_type,N_tot>result;
		for(size_type i=0;i<N_tot;++i)
		{
		    if(types[i]==1 or types[i]==-1)
		    {
			result[i]=1;
		    }
		    else if(types[i]==2 or types[i]==-2)
		    {
			result[i]=2;
		    }
		    else
		    {
			result[i]=0;
		    }
		}
		return result;
	    }

	    /* Public constructors: */
	    /*----------------------*/
	    
	    /* Non-allocating constructor: all colour index ranges are supposed
	     * to be N_c: */

	    colour_flow_QCD(const vector<object_type*,N_tot>& tensors):base_type(tensors),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    if(tensors[i]->rank()==2)
		    {
			const_weight*=2;
		    }
		}
	    }

	    /* Standalone-mode constructor. Argument array entry 0 is interpreted as a
	     * singlet, 1 or -1 as fundamental representation particles, and 2 for adjoint
	     * representations: */

	    colour_flow_QCD(const vector<int,N_tot>& types):base_type(make_rank_vector(types)),const_weight(1)
	    {
		for(size_type i=0;i<N_tot;++i)
		{
		    if(types[i]==2 or types[i]==-2)
		    {
			const_weight*=2;
		    }
		}
	    }

	    /* Public modifiers: */
	    /*-------------------*/

	    /* Generator method implementation. Quarks are generated as in the
	     * adjoint_QCD class, whereas gluons are generated by generated
	     * uniformly a point on the real (N*N-1)-dimensional sphere and
	     * contracting this vector with the SU(N) fundamental-representation
	     * matrices. */

	    bool generate()
	    {
		size_type colsum_factor(1);
		for(size_type i=0;i<N_tot;++i)
		{
		    if(this->colour_rank(i)==1)
		    {
			q_gen.generate();
			for(size_type j=0;j<N_c;++j)
			{
			    this->colour_entry(i,j)=q_norm*std::complex<value_type>(qvec[j<<1],qvec[(j<<1)+1]);
			}
		    }
		    if(this->colour_rank(i)==2)
		    {
			colsum_factor*=2;
			g_gen.generate();
			size_type c=0;
			this->colour_entry(i,0)=std::complex<value_type>(0,0);
			for(size_type j=1;j<N_c;++j)
			{
			    for(size_type k=0;k<j;++k)
			    {
				this->colour_entry(i,k*N_c+j)=g_norm*std::complex<value_type>(gvec[c],-gvec[c+1])/(value_type)2;
				this->colour_entry(i,j*N_c+k)=g_norm*std::complex<value_type>(gvec[c],gvec[c+1])/(value_type)2;
				c+=2;
			    }
			    value_type factor=g_norm*gvec[c]/std::sqrt(value_type(2*j*(j+1)));
			    for(size_type k=0;k<j;++k)
			    {
				this->colour_entry(i,k*(N_c+1))+=factor;
			    }
			    this->colour_entry(i,j*(N_c+1))=-(value_type)j*factor;
			    ++c;
			}
		    }
		}
		this->weight()=(value_type)const_weight/(value_type)colsum_factor;
		return true;
	    }

	    /* Weight evaluation method. Assigns the constant weight: */
	    
	    bool evaluate_weight()
	    {
		this->weight()=(value_type)const_weight;
		return true;
	    }

	    /* Public const methods: */
	    /*-----------------------*/

	    /* Clone method. */

	    colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>* clone() const
	    {
		return new colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>(*this);
	    }

	    /* Serialization: */
	    /*----------------*/

	    /* Polymorphic type identifier: */
	    
	    std::string type() const
	    {
		return "qcd";
	    }

	private:

	    /* Private static data: */
	    /*----------------------*/

	    /* Utility vectors: */

	    static vector<value_t,2*N_c>qvec;
	    static vector<value_t,N_c*N_c-1>gvec;

	    /* Helper generators: */

	    static uniform_sphere<value_t,2*N_c-1,rng_t>q_gen;
	    static uniform_sphere<value_t,N_c*N_c-2,rng_t>g_gen;

	    /* Normalisation factors: */

	    static const value_t q_norm;
	    static const value_t g_norm;

	    /* Private data: */
	    /*---------------*/

	    /* Constant weight copy: */

	    size_type const_weight;
    };
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const std::size_t colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::N_tot;
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>vector<value_t,2*N_c> colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::qvec;
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>vector<value_t,N_c*N_c-1> colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::gvec;
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>uniform_sphere<value_t,2*N_c-1,rng_t> colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::q_gen(&colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::qvec);
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>uniform_sphere<value_t,N_c*N_c-2,rng_t> colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::g_gen(&colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::gvec);
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const value_t colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::q_norm=std::sqrt((value_t)N_c);
    template<class value_t,std::size_t N_in,std::size_t N_out,std::size_t N_c,class rng_t>const value_t colour_flow_QCD<value_t,N_in,N_out,N_c,rng_t,true>::g_norm=std::sqrt(value_t(N_c*N_c-1));
}

#endif /*CAMGEN_QCD_COLS_H_*/

