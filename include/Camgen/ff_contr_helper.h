//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_FF_CONTR_HELPER_H_
#define CAMGEN_FF_CONTR_HELPER_H_

#include <Camgen/eval.h>
#include <Camgen/group.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Helper class for the contracted structure constants colour structure. It  *
 * lists the nonzero components of a recursive relation to optimise the      *
 * contraction procedure.                                                    *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class ff_contr_helper
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

			value_type f;
			for(size_type a=0;a<group_t::rank;++a)
			{
			    for(size_type b=0;b<group_t::rank;++b)
			    {
				for(size_type c=0;c<group_t::rank;++c)
				{
				    for(size_type d=0;d<group_t::rank;++d)
				    {
					f=value_type(0,0);
					for(size_type e=0;e<group_t::rank;++e)
					{
					    f+=((group_type::structure_constant(e,a,b))*(group_type::structure_constant(e,c,d)));
					}
					if(f!=value_type(0,0))
					{
					    values.push_back(f);
					    jumps[0].push_back(a*Feynrule_t::sizes[I]);
					    jumps[1].push_back(b*Feynrule_t::sizes[J]);
					    jumps[2].push_back(c*Feynrule_t::sizes[K]);
					    jumps[3].push_back(d*Feynrule_t::sizes[L]);
					}
				    }
				}
			    }
			}

			/* Secondly, we subtract each jump entry by the previous one to
			 * obtain the offsets: */

			for(size_type i=0;i<4;++i)
			{
			    backsets[i]=jumps[i].back();
			    for(int j=size()-1;j!=0;--j)
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

		/* Nonzero colour structure values: */

		static const value_type& value(size_type i)
		{
		    CAMGEN_ERROR_IF((i>=values.size()),"requested value out of range");
		    return values[i];
		}

		/* Offsets to nonzero colour structure components: */

		static int jump(size_type i,size_type j)
		{
		    CAMGEN_ERROR_IF((i>=4),"requested jump index out of range");
		    CAMGEN_ERROR_IF((j>=jumps[i].size()),"requested jump number out of range");
		    return jumps[i][j];
		}

		/* First nonzero colour structure component: */

		static size_type reset(size_type i)
		{
		    CAMGEN_ERROR_IF((i>=4),"requested index out of range");
		    return backsets[i];
		}

		/* Output method: */

		static void print()
		{
		    initialise();
		    for(size_type i=0;i<values.size();++i)
		    {
			std::cout<<jumps[0][i]<<","<<jumps[1][i]<<","<<jumps[2][i]<<","<<jumps[3][i]<<"\t\t"<<values[i]<<std::endl;
		    }
		}
	    
	    private:

		/* Nonzero colour structure values: */

		static std::vector<value_type> values;
		
		/* Offsets to nonzero colour structure entries: */
		
		static std::vector<int>jumps[4];

		/* Offsets to first nonzero colour structure entries: */

		static size_type backsets[4];
		
		/* Initilisation tag: */
		
		static bool initialised;
	};
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>bool ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::initialised=false;
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>std::vector<typename ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::value_type> ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::values;
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>std::vector<int> ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::jumps[4];
	template<class group_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>typename ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::size_type ff_contr_helper<group_t,Feynrule_t,I,J,K,L>::backsets[4]={0,0,0,0};
    }
}

#endif /*CAMGEN_FF_CONTR_HELPER_H_*/

