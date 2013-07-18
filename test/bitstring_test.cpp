//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/license_print.h>
#include <Camgen/bspart.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Program for checking the algorithmic functions on bit strings and the *
 * partition of bit strings.                                             *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */


namespace Camgen
{
    namespace test_utils
    {
	template<std::size_t N>bool check_partition(unsigned m,unsigned n,std::ostream& os=std::cerr)
	{
	    typedef typename bit_string_partition<N>::partition partition;
	    typedef typename bit_string_partition<N>::partition_iterator partition_iterator;

	    /* If upper-right-triangle partitions are requested to be checked,
	     * return true: */

	    if(n>m)
	    {
		return true;
	    }

	    /* Initialisaing and sorting the partition table: */

	    bit_string_partition<N>::initialise(n);
	    bit_string_partition<N>::sort();

	    /* Copy the vector of partitions that need to be checked: */

	    partition p=bit_string_partition<N>::get_partition(m,n);
	    bit_string<N>b=bit_string_partition<N>::get_partition(m,0)[0][0];

	    /* Checking whether each partition is unique in the vector: */

	    partition_iterator it=std::adjacent_find(p.begin(),p.end());
	    if(it != p.end())
	    {
		os<<"double counted partition found:"<<std::endl;
		bit_string_partition<N>::print_entry(os,m,n);
		return false;
	    }
	    for(it=p.begin();it!=p.end();++it)
	    {
		bit_string<N>b2;

		/* Checking whether the strings add up to the required string: */

		for(unsigned i=0;i<it->size();++i)
		{
		    b2^=((*it)[i]);
		}
		if(b2 != b)
		{
		    os<<"partition ";
		    bit_string_partition<N>::print_entry(os,m,n);
		    os<<" does not add up to "<<b<<std::endl;
		    return false;
		}
		b2.reset();

		/* Checking whether the strings are mutually disjoint: */

		for(unsigned i=0;i<it->size();++i)
		{
		    for(unsigned j=i+1;j<it->size();++j)
		    {
			if(((*it)[i]&(*it)[j]).any())
			{
			    os<<"in partition "<<std::endl;
			    bit_string_partition<N>::print_entry(os,m,n);
			    os<<"bitstrings "<<(*it)[i]<<" and "<<(*it)[j]<<" are not mutually disjoint."<<std::endl;
			    return false;
			}
		    }
		}
	    }
	    return true;
	}

	template<int N>bool check_bitstring_constructor(std::ostream& os=std::cerr)
	{
	    os<<"Checking constructor for length "<<N<<"..........";
	    os.flush();
	    bit_string<N>check;
	    for(unsigned i=0;i<(1<<N)+2;++i)
	    {
		bit_string<N>b(i);
		if(b!=check)
		{
		    os<<"Error in constructor of bitstring "<<i<<":";
		    os<<b<<" not equal to "<<check<<std::endl;
		    return false;
		}
		++check;
	    }
	    os<<"..........done."<<std::endl;
	    return true;
	}

	template<int N>bool check_bitstring_conversion(std::ostream& os=std::cerr)
	{
	    os<<"Checking conversion to integer for length "<<N<<"..........";
	    os.flush();
	    bit_string<N>check;
	    for(unsigned i=0;i<(1<<N)+2;++i)
	    {
		if(check.to_integer()!=(i%(1<<N)))
		{
		    os<<"Error in conversion to integer: ";
		    os<<check.to_integer()<<" not equal to "<<i%(1<<N)<<std::endl;
		    return false;
		}
		++check;
	    }
	    os<<"..........done."<<std::endl;
	    return true;
	}

	template<int N>bool check_bitstring_comparison(std::ostream& os=std::cerr)
	{
	    os<<"Checking comparison operators for length "<<N<<"..........";
	    os.flush();
	    bit_string<N>check;
	    bit_string<N>b;
	    for(unsigned i=0;i<(1<<N);++i)
	    {
		b.reset();
		for(unsigned j=0;j<i;++j)
		{
		    if(!(b<check))
		    {
			os<<"Error in comparison operator: ";
			os<<b<<" not smaller than "<<check<<std::endl;
			return false;
		    }
		    ++b;
		}
		++check;
	    }
	    os<<"..........done."<<std::endl;
	    return true;
	}

	template<int N>bool check_partition_table(std::ostream& os=std::cerr)
	{
	    os<<"Checking bitstring partition table for length "<<N<<"...........";
	    os.flush();
	    bit_string_partition<N>::initialise();
	    for(unsigned n=0;n<N;++n)
	    {
		for(unsigned m=0;m<=n;++m)
		{
		    if(!check_partition<N>(n,m,os))
		    {
			return false;
		    }
		}
	    }
	    os<<"..........done."<<std::endl;
	    return true;
	}

	template<int N>bool check_bitstrings(std::ostream& os=std::cerr)
	{
	    if(!test_utils::check_bitstring_constructor<N>(std::cerr))
	    {
		return false;
	    }
	    if(!test_utils::check_bitstring_conversion<N>(std::cerr))
	    {
		return false;
	    }
	    if(!test_utils::check_bitstring_comparison<N>(std::cerr))
	    {
		return false;
	    }
	    if(!test_utils::check_partition_table<N>(std::cerr))
	    {
		return false;
	    }
	    return true;
	}
    }
}

using namespace Camgen;

int main()
{
    license_print::initialise();
    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing bit string arithmetics............................"<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    if(!test_utils::check_bitstrings<3>(std::cerr))
    {
	return 1;
    }
    if(!test_utils::check_bitstrings<4>(std::cerr))
    {
	return 1;
    }
    if(!test_utils::check_bitstrings<5>(std::cerr))
    {
	return 1;
    }
    if(!test_utils::check_bitstrings<6>(std::cerr))
    {
	return 1;
    }
    if(!test_utils::check_bitstrings<7>(std::cerr))
    {
	return 1;
    }
    if(!test_utils::check_bitstrings<8>(std::cerr))
    {
	return 1;
    }
    if(!test_utils::check_bitstrings<9>(std::cerr))
    {
	return 1;
    }
    if(!test_utils::check_bitstrings<10>(std::cerr))
    {
	return 1;
    }
}

