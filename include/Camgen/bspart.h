//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_BSPART_H_
#define CAMGEN_BSPART_H_

#include <vector>
#include <algorithm>
#include <Camgen/debug.h>
#include <Camgen/bit_string.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the bitstring partition class. The class      *
 * template has a static initialisation method which takes an argument integer *
 * r, denoting the maximum number of bitstrings to participate in a partition. *
 * Subsequently it lists all (ordered) partitions of the bitstrings B_i, i<=r, *
 * where B_i denotes the string with the first i bits set and the rest zero.   *
 * After initialisation, any bit string can be quickly partitioned by          *
 * convoluting with the appropriate partitions.                                *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<std::size_t N>class bit_string_partition
    {
	public:
	    
	    /* Some useful type definitions: */

	    typedef std::size_t size_type;
	    typedef std::vector< std::vector< bit_string<N> > > partition;
	    typedef typename partition::iterator partition_iterator;
	    typedef typename std::vector< bit_string<N> > bit_string_iterator;
	    typedef typename std::vector< bit_string<N> > reverse_bit_string_iterator;

	    /* Initialisation phase: */

	    static void initialise(size_type lvl=N)
	    {
		/* Proceed with 'N' for unphysical levels: */
		
		if(lvl>N)
		{
		    lvl=N;
		}
		
		/* If the required level is larger than the already generated
		 * partition level, proceed: */
		
		if(lvl>max_lvl)
		{
		    /* If it is the first nontrivial initialisation, fill the
		     * first column with the trivial partitions: */
		    
		    if(max_lvl==0)
		    {
			bit_string<N>b;
			for(size_type i=0;i<N;++i)
			{
			    b.set(i);
			    partition_table[i][0].resize(1);
			    partition_table[i][0][0].resize(1);
			    partition_table[i][0][0][0]=b;
			}
		    }

		    /* If the required splitting level is larger than one,
		     * proceed: */
		    
		    if(lvl>1)
		    {
			/* If the [1][1]-partition was not yet generated, do it
			 * manually: */

			if(max_lvl<2)
			{
			    partition_table[1][1].resize(1);
			    partition_table[1][1][0].resize(2);
			    partition_table[1][1][0][0].set(0);
			    partition_table[1][1][0][1].set(1);
			}
			for(size_type i=2;i<N;++i)
			{
			    /* Recursively construct the other partitions: */

			    for(size_type j=max_lvl;j<i;++j)
			    {
				/* Fill the lower triangle of the partition
				 * table, starting from the present level, and
				 * ending at the required level: */

				if(j<lvl)
				{
				    /* Resize the table entries by the recursion
				     * relation obeyed by the Stirling numbers:
				     * */

				    partition_table[i][j].resize(partition_table[i-1][j-1].size()+(j+1)*partition_table[i-1][j].size());
				    
				    /* Construct the partition from the ones in
				     * the previous row: */

				    size_type n=0;
				    for(size_type l=0;l<partition_table[i-1][j].size();++l)
				    {
					for(int m=j;m!=-1;--m)
					{
					    partition_table[i][j][n]=partition_table[i-1][j][l];
					    partition_table[i][j][n][m].set(i);
					    ++n;
					}
				    }
				    for(size_type l=0;l<partition_table[i-1][j-1].size();++l)
				    {
					partition_table[i][j][n]=partition_table[i-1][j-1][l];
					partition_table[i][j][n].resize(partition_table[i][j][n].size()+1);
					partition_table[i][j][n].back().set(i);
					++n;
				    }
				}
				else
				{
				    break;
				}
			    }
			    if(i<lvl)
			    {
				/* Fill the diagonal entries by the maximal
				 * partition: */

				partition_table[i][i].resize(1);
				partition_table[i][i][0].resize(i+1);
				for(size_type j=0;j<i+1;++j)
				{
				    partition_table[i][i][0][j].set(j);
				}
			    }
			}
		    }

		    /* Copy the level to the maximum level: */

		    max_lvl=lvl;
		}
	    }

	    /* Static function sorting all partitions: */

	    static void sort()
	    {
		if(!sorted)
		{
		    for(size_type i=0;i<N;++i)
		    {
			for(size_type j=0;j<max_lvl;++j)
			{
			    for(size_type k=0;k<partition_table[i][j].size();++k)
			    {
				std::sort(partition_table[i][j][k].begin(),partition_table[i][j][k].end());
			    }
			    std::sort(partition_table[i][j].begin(),partition_table[i][j].end());
			}
		    }
		    sorted=true;
		}
	    }

	    /* Print the number of partitions in all lower-triangle table
	     * entries: */

	    static std::ostream& print_sizes(std::ostream& os)
	    {
		for(size_type i=0;i<N;++i)
		{
		    for(size_type j=0;j<N;++j)
		    {
			os.width(6);
			os<<partition_table[i][j].size()<<"   ";
		    }
		    os<<std::endl;
		}
		return os;
	    }

	    /* Print a particular partition: */

	    static std::ostream& print_entry(std::ostream& os,size_type n,size_type k)
	    {
		os<<"[";
		if(n<N and k<N)
		{
		    size_type s=partition_table[n][k].size();
		    if(s>0)
		    {
			for(size_type i=0;i<s-1;++i)
			{
			    for(size_type j=0;j<partition_table[n][k][i].size();++j)
			    {
				os<<"["<<partition_table[n][k][i][j]<<"]";
			    }
			    os<<",";
			}
			for(size_type j=0;j<partition_table[n][k][s-1].size();++j)
			{
			    os<<"["<<partition_table[n][k][s-1][j]<<"]";
			}
		    }
		}
		os<<"]";
		return os;
	    }

	    /* Construct all partitions of a bitstring b into n parts, by
	     * convoluting with the appropriate table entries: */

	    static partition make_partition(const bit_string<N>& b,size_type n)
	    {
		partition result;
		size_type c=b.count();
		if(n>c)
		{
		    return result;
		}
		if(n>max_lvl)
		{
		    return result;
		}
		result.resize(partition_table[c-1][n-1].size());
		for(size_type i=0;i<result.size();++i)
		{
		    result[i].resize(n,b);
		    for(size_type j=0;j<n;++j)
		    {
			result[i][j]*=(partition_table[c-1][n-1][i][j]);
		    }
		}
		return result;
	    }

	    /* Pass partition entry: */

	    static const partition& get_partition(size_type i,size_type j)
	    {
		CAMGEN_ERROR_IF((i>=N),"First partition entry out of range");
		CAMGEN_ERROR_IF((j>=N),"Second partition entry out of range");
		return partition_table[i][j];
	    }

	    /* Pass maximum partition level: */

	    static size_type max_level()
	    {
		return max_lvl;
	    }

	    /* Pass sorted boolean: */

	    static bool is_sorted()
	    {
		return sorted;
	    }

	private:

	    /* Partition table static data member. The rows specify the number
	     * of bits set in the partitioned string, the columns the number of
	     * strings in a partition. Hence the upper-right triangle in the
	     * table consists of empty partitions: */

	    static partition partition_table[N][N];
	    
	    /* Current maximal partitioned level: */
	    
	    static size_type max_lvl;

	    /* Boolean denoting whether the partitions are sorted: */

	    static bool sorted;
    };
    template<std::size_t N>typename bit_string_partition<N>::partition bit_string_partition<N>::partition_table[N][N];
    template<std::size_t N>typename bit_string_partition<N>::size_type bit_string_partition<N>::max_lvl=0;
    template<std::size_t N>bool bit_string_partition<N>::sorted=false;
}


#endif /*CAMGEN_BSPART_H_*/


