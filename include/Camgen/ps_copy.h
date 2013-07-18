//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file ps_copy.h
    \brief Classes to copy phase space variables between different CM_algorithm instances.
 */

#ifndef CAMGEN_PS_COPY_H_
#define CAMGEN_PS_COPY_H_

#include <Camgen/CM_algo.h>

namespace Camgen
{
    template<class model_t1,class model_t2>class phase_space_copy;

    /* Momentum copy helper class, issues a warning by default: */

    template<class model_t1,class model_t2,std::size_t D1=model_t1::dimension,std::size_t D2=model_t2::dimension>class momentum_copy
    {
	private:

	    friend class phase_space_copy<model_t1,model_t2>;
	    
	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"no implementation defined for unequal-sized vector copy"<<endlog;
	    }
    };

    /* Momentum copy helper class template specialisation for same-dimensional
     * models: */

    template<class model_t1,class model_t2,std::size_t D>class momentum_copy<model_t1,model_t2,D,D>
    {
	private:
	    
	    friend class phase_space_copy<model_t1,model_t2>;
	    
	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		if(D!=0)
		{
		    for(std::size_t i=0;i<N_in+N_out;++i)
		    {
			if(it1->get_phase_space(i)!=NULL and it2->get_phase_space(i)!=NULL)
			{
			    it2->get_phase_space(i)->momentum()=it1->get_phase_space(i)->momentum();
			}
		    }
		}
	    }
    };

    /* Helicity copy class template, issues a warning by default: */

    template<class model_t1,class model_t2,bool q1=model_t1::continuous_helicities,bool q2=model_t2::continuous_helicities>class helicity_copy
    {
	private:
	    
	    friend class phase_space_copy<model_t1,model_t2>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"no implementation defined for copying continuous to discrete helicity variables"<<endlog;
	    }
    };

    /* Helicity copy class template specialisation for discrete helicity types:
     * */

    template<class model_t1,class model_t2>class helicity_copy<model_t1,model_t2,false,false>
    {
	private:
	    
	    friend class phase_space_copy<model_t1,model_t2>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		for(std::size_t i=0;i<N_in+N_out;++i)
		{
		    if(it1->get_phase_space(i)!=NULL and it2->get_phase_space(i)!=NULL)
		    {
			if(it1->get_phase_space(i)->helicities()==it2->get_phase_space(i)->helicities())
			{
			    it2->get_phase_space(i)->helicity()=it1->get_phase_space(i)->helicity();
			}
			else
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particles of different spin ignored in helicity copy function"<<endlog;
			}
		    }
		}
	    }
    };

    /* Helicity copy class template specialisation for discrete-to-continuous helicity copying: */

    template<class model_t1,class model_t2>class helicity_copy<model_t1,model_t2,false,true>
    {
	private:
	    
	    friend class phase_space_copy<model_t1,model_t2>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		for(std::size_t i=0;i<N_in+N_out;++i)
		{
		    if(it1->get_phase_space(i)!=NULL and it2->get_phase_space(i)!=NULL)
		    {
			if(it1->get_phase_space(i)->helicities()==it2->get_phase_space(i)->helicities())
			{
			    for(int j=-it2->get_phase_space(i)->maximal_helicity();j<0;++j)
			    {
				it2->get_phase_space(i)->helicity_phase(j)=std::complex<typename model_t2::value_type>(0,0);
			    }
			    if(it2->get_phase_space(i)->has_zero_helicity())
			    {
				it2->get_phase_space(i)->helicity_phase(0)=std::complex<typename model_t2::value_type>(0,0);
			    }
			    for(int j=1;j<=it2->get_phase_space(i)->maximal_helicity();++j)
			    {
				it2->get_phase_space(i)->helicity_phase(j)=std::complex<typename model_t2::value_type>(0,0);
			    }
			    it2->get_phase_space(i)->helicity_phase(it1->get_phase_space(i)->helicity()).real()=(typename model_t2::value_type)0;
			}
			else
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particles of different spin ignored in helicity copy function"<<endlog;
			}
		    }
		}
	    }
    };

    /* Helicity copy class template specialisation for continuous helicity types:
     * */

    template<class model_t1,class model_t2>class helicity_copy<model_t1,model_t2,true,true>
    {
	private:
	    
	    friend class phase_space_copy<model_t1,model_t2>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		for(std::size_t i=0;i<N_in+N_out;++i)
		{
		    if(it1->get_phase_space(i)!=NULL and it2->get_phase_space(i)!=NULL)
		    {
			if(it1->get_phase_space(i)->helicities()==it2->get_phase_space(i)->helicities())
			{
			    for(int j=-it2->get_phase_space(i)->maximal_helicity();j<0;++j)
			    {
				it2->get_phase_space(i)->helicity_phase(j)=it1->get_phase_space(i)->helicity_phase(j);
			    }
			    if(it2->get_phase_space(i)->has_zero_helicity())
			    {
				it2->get_phase_space(i)->helicity_phase(0)=it1->get_phase_space(i)->helicity_phase(0);
			    }
			    for(int j=1;j<=it2->get_phase_space(i)->maximal_helicity();++j)
			    {
				it2->get_phase_space(i)->helicity_phase(j)=it1->get_phase_space(i)->helicity_phase(j);
			    }
			}
			else
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"particles of different spin ignored in helicity copy function"<<endlog;
			}
		    }
		}
	    }
    };

    /* Colour copy class template helper: */

    template<class model_t1,class model_t2,bool q=(model_t1::coloured and model_t2::coloured)>class colour_copy;
    
    /* Colour copy class helper helper class: */

    template<class model_t1,class model_t2,bool q1=model_t1::continuous_colours,bool q2=model_t2::continuous_colours>class colour_copy_helper
    {
	private:

	    friend class colour_copy<model_t1,model_t2,true>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		log(log_level::warning)<<CAMGEN_STREAMLOC<<"no implementation defined for copying continuous to discrete helicity variables"<<endlog;
	    }
    };

    /* Colour copy helper helper class template specialisation for discrete
     * colour models: */

    template<class model_t1,class model_t2>class colour_copy_helper<model_t1,model_t2,false,false>
    {
	private:

	    friend class colour_copy<model_t1,model_t2,true>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		for(std::size_t i=0;i<N_in+N_out;++i)
		{
		    if(it1->get_phase_space(i)!=NULL and it2->get_phase_space(i)!=NULL)
		    {
			if(it1->get_phase_space(i)->particle_type->get_colour_index_ranges()==it2->get_phase_space(i)->particle_type->get_colour_index_ranges())
			{
			    for(std::size_t a=0;a<it2->get_phase_space(i)->colour_rank();++a)
			    {
				if(it1->get_phase_space(i)->colour_range(a)==it2->get_phase_space(i)->colour_range(a))
				{
				    it2->get_phase_space(i)->colour(a)=it1->get_phase_space(i)->colour(a);
				}
				else
				{
				    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to copy unequal-range colour degrees of ignored"<<endlog;
				}
			    }
			}
			else
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to copy unequal-rank colour degrees of freedom ignored"<<endlog;
			}
		    }
		}
	    }
    };

    /* Colour copy helper helper class template specialisation for continuous
     * colour models: */

    template<class model_t1,class model_t2>class colour_copy_helper<model_t1,model_t2,true,true>
    {
	private:

	    friend class colour_copy<model_t1,model_t2,true>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		for(std::size_t i=0;i<N_in+N_out;++i)
		{
		    if(it1->get_phase_space(i)!=NULL and it2->get_phase_space(i)!=NULL)
		    {
			if(it1->get_phase_space(i)->particle_type->get_colour_index_ranges()==it2->get_phase_space(i)->particle_type->get_colour_index_ranges())
			{
			    for(std::size_t a=0;a<it2->get_phase_space(i)->get_colours()->size();++a)
			    {
				it2->get_phase_space(i)->colour_entry(a)=it1->get_phase_space(i)->colour_entry(a);
			    }
			}
			else
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"attempt to copy unequal-size colour tensors ignored"<<endlog;
			}
		    }
		}
	    }
    };

    /* Colour copy helper helper class template specialisation for discrete-to-continuous
     * colour copying: */

    template<class model_t1,class model_t2>class colour_copy_helper<model_t1,model_t2,false,true>
    {
	private:

	    friend class colour_copy<model_t1,model_t2,true>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		typedef typename model_t2::value_type value_type;
		for(std::size_t i=0;i<N_in+N_out;++i)
		{
		    if(it1->get_phase_space(i)!=NULL and it2->get_phase_space(i)!=NULL)
		    {
			if(it1->get_phase_space(i)->particle_type->get_colour_index_ranges()==it2->get_phase_space(i)->particle_type->get_colour_index_ranges())
			{
			    it2->get_phase_space(i)->get_colours()->reset();
			    (*(it2->get_phase_space(i)->get_colours()))(*(it1->get_colours()))=std::complex<value_type>(1,0);
			}
			else
			{
			    log(log_level::warning)<<CAMGEN_STREAMLOC<<"Attempt to copy unequal-size colour tensors ignored"<<endlog;
			}
		    }
		}
	    }
    };

    /* Colour copy class template specialisation for colourless models: */

    template<class model_t1,class model_t2>class colour_copy<model_t1,model_t2,false>
    {
	private:
	    
	    friend class phase_space_copy<model_t1,model_t2>;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2){}
    };

    /* Colour copy class template specialisation for colourful models: */

    template<class model_t1,class model_t2>class colour_copy<model_t1,model_t2,true>
    {
	private:
	    
	    friend class phase_space_copy<model_t1,model_t2>;
	    
	    typedef colour_copy_helper<model_t1,model_t2,model_t1::continuous_colours,model_t2::continuous_colours> helper_class;

	    template<std::size_t N_in,std::size_t N_out>static void apply(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		helper_class::template apply<N_in,N_out>(it1,it2);
	    }
    };

    /* Phase space copy wrapper class: */

    template<class model_t1,class model_t2>class phase_space_copy
    {
	public:

	    typedef model_t1 model_type1;
	    typedef model_t2 model_type2;
	    typedef momentum_copy<model_t1,model_t2,model_t1::dimension,model_t2::dimension> momentum_cp_type;
	    typedef helicity_copy<model_t1,model_t2,model_t1::continuous_helicities,model_t2::continuous_helicities> helicity_cp_type;
	    typedef colour_copy<model_t1,model_t2,model_t1::coloured and model_t2::coloured> colour_cp_type;

	    template<std::size_t N_in,std::size_t N_out>static void momenta(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		momentum_cp_type::template apply<N_in,N_out>(it1,it2);
	    }

	    template<std::size_t N_in,std::size_t N_out>static void helicities(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		helicity_cp_type::template apply<N_in,N_out>(it1,it2);
	    }

	    template<std::size_t N_in,std::size_t N_out>static void colours(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
	    {
		colour_cp_type::template apply<N_in,N_out>(it1,it2);
	    }
    };

    /// Function copying momenta from CM tree iterator it1 to it2.

    template<class model_t1,class model_t2,std::size_t N_in,std::size_t N_out>void copy_momenta(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
    {
	phase_space_copy<model_t1,model_t2>::template momenta<N_in,N_out>(it1,it2);
    }

    /// Function copying helicities from CM tree iterator it1 to it2.

    template<class model_t1,class model_t2,std::size_t N_in,std::size_t N_out>void copy_helicities(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
    {
	phase_space_copy<model_t1,model_t2>::template helicities<N_in,N_out>(it1,it2);
    }

    /// Function copying colours from CM tree iterator it1 to it2.

    template<class model_t1,class model_t2,std::size_t N_in,std::size_t N_out>void copy_colours(typename CM_algorithm<model_t1,N_in,N_out>::const_tree_iterator it1,typename CM_algorithm<model_t2,N_in,N_out>::tree_iterator it2)
    {
	phase_space_copy<model_t1,model_t2>::template colours<N_in,N_out>(it1,it2);
    }
}

#endif /*CAMGEN_PS_COPY_H_*/

