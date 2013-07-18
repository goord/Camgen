//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_TT_CONTR_HELPER_H_
#define CAMGEN_TT_CONTR_HELPER_H_

#include <Camgen/eval.h>
#include <Camgen/rep.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Helper class for the TT_contr-vertex colour structure. It lists the nonzero *
 * group-index contracted generator product components in a recursion relation *
 * to optimise the contraction procedure.                                      *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>class TT_contr_helper
	{
	    public:

		/* Useful type definitions: */

		typedef typename evaluate<Feynrule_t>::size_type size_type;
		typedef typename evaluate<Feynrule_t>::value_type value_type;
		
		/* Representation type, providing the generator matrices: */
		
		typedef typename rep_t::template implementation<typename evaluate<Feynrule_t>::model_type::value_type> rep_type;

		/* Initialisation phase, creating the vector of 'jumps' of the iterators
		 * to nonzero entries: */

		static void initialise()
		{
		    if(!initialised)
		    {
			size_type Isize=Feynrule_t::sizes[I];
			size_type Jsize=Feynrule_t::sizes[J];
			size_type Ksize=Feynrule_t::sizes[K];
			size_type Lsize=Feynrule_t::sizes[L];

			rep_type::initialise();

			/* First, list the nonzero components in the 'jumps' vectors: */

			value_type temp;
			for(size_type i=0;i<rep_t::dimension;++i)
			{
			    for(size_type j=0;j<rep_t::dimension;++j)
			    {
				for(size_type k=0;k<rep_t::dimension;++k)
				{
				    for(size_type l=0;l<rep_t::dimension;++l)
				    {
					temp=value_type(0,0);
					for(size_type a=0;a<rep_t::group_type::rank;++a)
					{
					    temp+=((rep_type::generator(a,i,j))*(rep_type::generator(a,k,l)));
					}
					if(temp != value_type(0,0))
					{
					    values.push_back(temp);
					    jumps[0].push_back(i*Isize);
					    jumps[1].push_back(j*Jsize);
					    jumps[2].push_back(k*Ksize);
					    jumps[3].push_back(l*Lsize);
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
		    CAMGEN_ERROR_IF((i>=values.size()),"requested value index out of range");
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
			std::cout<<jumps[0][i]<<","<<jumps[1][i]<<","<<jumps[1][i]<<","<<jumps[3][i]<<"\t\t"<<values[i]<<std::endl;
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
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>bool TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::initialised=false;
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>std::vector<typename TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::value_type> TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::values;
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>std::vector<int> TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::jumps[4];
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K,std::size_t L>typename TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::size_type TT_contr_helper<rep_t,Feynrule_t,I,J,K,L>::backsets[4]={0,0,0,0};
    }
}

#endif /*CAMGEN_TT_CONTR_HELPER_H_*/

