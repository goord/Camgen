//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <iostream>
#include <Camgen/license_print.h>
#include <Camgen/num_utils.h>
#include <Camgen/vector.h>
#include <Camgen/Lorentz.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Euclidean.h>
#include <Camgen/stdrand.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Program for checking the basic Lorentz vector manipulation functions: *
 * boost from rest frame, boost to rest frame, etc.                      *
 *                                                                       *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /* Helper class for checking spacetime implementation: */

    template<class value_t, std::size_t D,class spacetime_t>class spacetime_test_helper
    {
	public:

	    typedef vector<value_t,D> vector_type;
	    typedef typename vector_type::size_type size_type;
	    typedef typename vector_type::value_type value_type;
	    typedef typename spacetime_t::template implementation<value_t,D> spacetime_type;
	    typedef typename spacetime_type::base_type spacetime_base_type;

	    /* Checks symmetry and diagonality of the metric: */

	    static bool check_metric()
	    {
		spacetime_type::initialise();

		for(size_type i=0;i<D;++i)
		{
		    for(size_type j=0;j<i;++j)
		    {
			// check if metric is symmetric:
			if(spacetime_type::metric(i,j)!=spacetime_type::metric(j,i))
			{
			    return false;
			}
			// check if off-diagonal elements vanish if declared diagonal:
			if(spacetime_type::diagonal and spacetime_type::metric(i,j)!=(value_type)0)
			{
			    return false;
			}
		    }
		}
		return true;
	    }

	    /* Checks symmetry and consistency of the dot product with the metric: */

	    static bool check_dot_product(const vector_type& u,const vector_type& v)
	    {
		// check whether dot product is (exactly) symmetric:
		if(spacetime_type::dot(u,v)!=spacetime_type::dot(v,u))
		{
		    return false;
		}
		// check whether dot product is correctly overloaded:
		if(!equals(spacetime_base_type::dot(u,v),spacetime_type::dot(u,v)))
		{
		    return false;
		}
		return true;
	    }

	    static bool check_dot_product_loop(size_type N,const value_type& scale)
	    {
		vector_type u,v;
		for(size_type i=0;i<N;++i)
		{
		    for(size_type mu=0;mu<D;++mu)
		    {
			u[mu]=random_number_stream<value_type,std::random>::throw_number(-scale,scale);
			v[mu]=random_number_stream<value_type,std::random>::throw_number(-scale,scale);
		    }
		    if(!check_dot_product(u,v))
		    {
			return false;
		    }
		}
		return true;
	    }
    };

    /* Helper class checking the Lorentz transformations: */

    template<class value_t, std::size_t D>class Lorentz_test_helper
    {
	public:

	    typedef vector<value_t,D> vector_type;
	    typedef typename vector_type::size_type size_type;
	    typedef typename vector_type::value_type value_type;
	    typedef vector<value_t,D-1> space_vector_type;
	    typedef typename Minkowski_type::template implementation<value_t,D> spacetime_type;

	    static bool check_Lorentz_transform_loop(size_type m, size_type n,const value_type& scale)
	    {
		vector_type u,v;
		value_type sv=0;
		for(size_type i=0;i<m;++i)
		{
		    do
		    {
			for(size_type mu=0;mu<D;++mu)
			{
			    v[mu]=random_number_stream<value_type,std::random>::throw_number(-scale,scale);
			}
			sv=spacetime_type::dot(v,v);
		    }
		    while(sv<=(value_type)0);

		    if(!check_Lorentz_transform(v))
		    {
			return false;
		    }
		    for(size_type j=0;j<n;++j)
		    {
			for(size_type mu=0;mu<D;++mu)
			{
			    u[mu]=random_number_stream<value_type,std::random>::throw_number(-scale,scale);
			}
			if(!check_Lorentz_transform(u,v))
			{
			    return false;
			}
		    }
		}
		return true;
	    }

	    static bool check_Lorentz_transform(const vector_type& u,const vector_type& v)
	    {
		vector_type w=u;
		boost_to_restframe(w,v);
		vector_type w2=copy_boost_to_restframe(u,v);
		
		// check whether copy_boost correctly copies:
		if(!equal_sequences(w,w2))
		{
		    return false;
		}

		// check whether invariant mass is the same:
		if(!equals(spacetime_type::dot(w,w),spacetime_type::dot(u,u)))
		{
		    return false;
		}

		vector_type w3=copy_boost_from_restframe(w,v);

		// check whether invariant mass is the same:
		if(!equals(spacetime_type::dot(w3,w3),spacetime_type::dot(w,w)))
		{
		    return false;
		}
		
		// check whether boost_from is the inverse of boost_to:
		if(!equal_sequences(w3,u))
		{
		    return false;
		}

		return true;
	    }

	    static bool check_Lorentz_transform(const vector_type& v)
	    {
		vector_type v2=v;
		boost_to_restframe(v2,v);

		// check whether invariant mass is the same:
		if(!equals(spacetime_type::dot(v2,v2),spacetime_type::dot(v,v)))
		{
		    return false;
		}

		space_vector_type v2vec=make_space_vector(v2);

		//check whether boost to self-restframe gives zero 3-vector:
		for(size_type i=0;i<v2vec.size();++i)
		{
		    if(!equals(v2vec[i],(value_type)0))
		    {
			return false;
		    }
		}

		vector_type v3=v2;
		boost_from_restframe(v3,v);
		
		// check whether boost_from is the inverse of boost_to:
		if(!equal_sequences(v3,v))
		{
		    return false;
		}

		return true;
	    }
	    

	    /* Returns the 3-vector correspoding to the argument Lorentz vector: */

	    static space_vector_type make_space_vector(const vector_type& v)
	    {
		space_vector_type w;

		for(size_type mu=0;mu<spacetime_type::timelike_direction;mu++)
		{
		    w[mu]=v[mu];
		}
		for(size_type mu=spacetime_type::timelike_direction;mu<D-1;mu++)
		{
		    w[mu]=v[mu+1];
		}
		return w;
	    }
    };
}


using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing Lorentz vectors..................................."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;

    std::cout<<"Checking whether 4-dimensional Minkowski space has consistent metric..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,4,Minkowski_type>::check_metric())
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 4-dimensional Minkowski space has consistent dot product..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,4,Minkowski_type>::check_dot_product_loop(1000000,1e+6))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 10-dimensional Minkowski space has consistent metric..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,10,Minkowski_type>::check_metric())
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 10-dimensional Minkowski space has consistent dot product..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,10,Minkowski_type>::check_dot_product_loop(1000000,1e+6))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 3-dimensional Euclidean space has consistent metric..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,3,Euclidean_type>::check_metric())
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 3-dimensional Euclidean space has consistent dot product..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,3,Euclidean_type>::check_dot_product_loop(1000000,1e+6))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
    
    std::cout<<"Checking whether 10-dimensional Euclidean space has consistent metric..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,10,Euclidean_type>::check_metric())
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether 10-dimensional Euclidean space has consistent dot product..........";
    std::cout.flush();
    if(!spacetime_test_helper<double,10,Euclidean_type>::check_dot_product_loop(1000000,1e+6))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether Lorentz transform of 4-vectors is consistent..........";
    std::cout.flush();
    if(!Lorentz_test_helper<double,4>::check_Lorentz_transform_loop(10000,100,1e+6))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether Lorentz transform of 7-vectors is consistent..........";
    std::cout.flush();
    if(!Lorentz_test_helper<double,7>::check_Lorentz_transform_loop(10000,100,1e+6))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;

    std::cout<<"Checking whether Lorentz transform of 10-vectors is consistent..........";
    std::cout.flush();
    if(!Lorentz_test_helper<double,10>::check_Lorentz_transform_loop(10000,100,1e+6))
    {
	return 1;
    }
    std::cout<<"..........done."<<std::endl;
}

