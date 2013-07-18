//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_F_HELPER_H_
#define CAMGEN_F_HELPER_H_

#include <Camgen/eval.h>
#include <Camgen/group.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Helper class for the f-vertex colour structrure. It lists the nonzero     *
 * components of a recursive relation involving group structure constants to *
 * optimise the contraction procedure.                                       *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>class f_helper
	{
	    public:
		/* Useful type definitions: */

		typedef typename evaluate<Feynrule_t>::value_type value_type;
		typedef typename evaluate<Feynrule_t>::size_type size_type;
		
		/* Group type, the source of the structure constants: */
		
		typedef typename group_t::template implementation<typename evaluate<Feynrule_t>::model_type::value_type> group_type;

		/* Initialisation phase, creating the vector of 'jumps' of the iterators
		 * to nonzero entries: */

		static void initialise()
		{
		    if(!initialised)
		    {
			group_type::initialise();

			/* First, list the nonzero components in the 'jumps' vectors: */

			for(size_type a=0;a<group_t::rank;++a)
			{
			    for(size_type b=0;b<group_t::rank;++b)
			    {
				for(size_type c=0;c<group_t::rank;++c)
				{
				    if(group_type::structure_constant(a,b,c) != value_type(0,0))
				    {
					values.push_back(group_type::structure_constant(a,b,c));
					jumps[0].push_back(a*Feynrule_t::sizes[I]);
					jumps[1].push_back(b*Feynrule_t::sizes[J]);
					jumps[2].push_back(c*Feynrule_t::sizes[K]);
				    }
				}
			    }
			}

			/* Secondly, we subtract each jump entry by the previous one to
			 * obtain the offsets: */

			for(size_type i=0;i<3;++i)
			{
			    backsets[i]=jumps[i].back();
			    for(int j=size()-1;j != 0;--j)
			    {
				jumps[i][j]-=jumps[i][j-1];
			    }
			}
			initialised=true;
		    }
		}
		
		/* Number of nonzero components: */

		static size_type size()
		{
		    return values.size();
		}

		/* Nonzero structure constant values: */

		static const value_type& value(size_type i)
		{
		    CAMGEN_ERROR_IF((i>=values.size()),"requested value out of range");
		    return values[i];
		}

		/* Offsets to nonzero structure constants: */

		static int jump(size_type i,size_type j)
		{
		    CAMGEN_ERROR_IF((i>=3),"requested jump index out of range");
		    CAMGEN_ERROR_IF((j>=jumps[i].size()),"requested jump number out of range");
		    return jumps[i][j];
		}

		/* First nonzero structure constant index: */

		static size_type reset(size_type i)
		{
		    CAMGEN_ERROR_IF((i>=3),"requested index out of range");
		    return backsets[i];
		}

		/* Output method: */

		static void print()
		{
		    initialise();
		    for(size_type i=0;i<values.size();++i)
		    {
			std::cout<<jumps[0][i]<<","<<jumps[1][i]<<","<<jumps[2][i]<<"\t\t"<<values[i]<<std::endl;
		    }
		}
	    
	    private:

		/* Nonzero structure constant values: */

		static std::vector<value_type> values;
		
		/* Offsets to nonzero structure constants: */
		
		static std::vector<int>jumps[3];

		/* Offsets to first nonzero structure constants: */

		static size_type backsets[3];
		
		/* Initilisation tag: */
		
		static bool initialised;
	};
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>bool f_helper<group_t,Feynrule_t,I,J,K>::initialised=false;
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>std::vector<typename f_helper<group_t,Feynrule_t,I,J,K>::value_type> f_helper<group_t,Feynrule_t,I,J,K>::values;
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>std::vector<int> f_helper<group_t,Feynrule_t,I,J,K>::jumps[3];
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>typename f_helper<group_t,Feynrule_t,I,J,K>::size_type f_helper<group_t,Feynrule_t,I,J,K>::backsets[3]={0,0,0};
    }
}

#endif /*CAMGEN_F_HELPER_H_*/

