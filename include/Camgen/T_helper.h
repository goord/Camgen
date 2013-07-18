//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_T_HELPER_H_
#define CAMGEN_T_HELPER_H_

#include <Camgen/eval.h>
#include <Camgen/rep.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
 * Helper class for the T-vertex colour structure. It lists the nonzero      *
 * generator components in a recursive relation to optimise the contraction  *   
 * procedure                                                                 *
 *                                                                           *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace colour_tensor
    {
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>class T_helper
	{
	    public:

		/* Useful type definitions: */

		typedef typename evaluate<Feynrule_t>::value_type value_type;
		typedef typename evaluate<Feynrule_t>::size_type size_type;
		
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
			rep_type::initialise();

			/* First, list the nonzero components in the 'jumps' vectors: */

			for(size_type a=0;a<rep_t::group_type::rank;++a)
			{
			    for(size_type i=0;i<rep_t::dimension;++i)
			    {
				for(size_type j=0;j<rep_t::dimension;++j)
				{
				    if(rep_type::generator(a,i,j) != (value_type)0)
				    {
					values.push_back(rep_type::generator(a,i,j));
					jumps[0].push_back(a*Isize);
					jumps[1].push_back(i*Jsize);
					jumps[2].push_back(j*Ksize);
				    }
				}
			    }
			}

			/* Secondly, we subtract each jump entry by the previous one to
			 * obtain the offsets: */

			for(size_type i=0;i<3;++i)
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

		/* Nonzero generator values: */

		static const value_type& value(size_type i)
		{
		    CAMGEN_ERROR_IF((i>=values.size()),"requested value out of range");
		    return values[i];
		}

		/* Offsets to nonzero generator components: */

		static int jump(size_type i,size_type j)
		{
		    CAMGEN_ERROR_IF((i>=3),"requested jump index out of range");
		    CAMGEN_ERROR_IF((j>=jumps[i].size()),"requested jump number out of range");
		    return jumps[i][j];
		}

		/* First nonzero generator component: */

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

		/* Nonzero colour structure values: */

		static std::vector<value_type> values;
		
		/* Offsets to nonzero colour structure entries: */
		
		static std::vector<int>jumps[3];

		/* Offsets to first nonzero colour structure entries: */

		static size_type backsets[3];
		
		/* Initilisation tag: */
		
		static bool initialised;
	};
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>bool T_helper<rep_t,Feynrule_t,I,J,K>::initialised=false;
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>std::vector<typename T_helper<rep_t,Feynrule_t,I,J,K>::value_type> T_helper<rep_t,Feynrule_t,I,J,K>::values;
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>std::vector<int> T_helper<rep_t,Feynrule_t,I,J,K>::jumps[3];
	template<class rep_t,class Feynrule_t,std::size_t I,std::size_t J,std::size_t K>typename T_helper<rep_t,Feynrule_t,I,J,K>::size_type T_helper<rep_t,Feynrule_t,I,J,K>::backsets[3]={0,0,0};
    }
}

#endif /*CAMGEN_T_HELPER_H_*/

