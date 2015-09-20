//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file uni_sphere.h
    \brief uniform multi-dimensional sphere generators.
 */


#ifndef CAMGEN_UNI_SPHERE_H_
#define CAMGEN_UNI_SPHERE_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Classes to generate uniformly points on the n-sphere. The code uses the *
 * Marsaglia algorithm, with optimized specialisations for the 1-,2- and 3-*
 * sphere generations.                                                     *
 *                                                                         *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/vector.h>
#include <Camgen/combs.h>
#include <Camgen/norm_gen.h>
#include <Camgen/obj_alloc.h>

namespace Camgen
{
    /// Uniform sphere generator class.

    template<class value_t,std::size_t D,class rng_t>class uniform_sphere: public object_allocator< vector<value_t,D+1> >, public MC_generator<value_t>
    {
	typedef object_allocator< vector<value_t,D+1> > base_type;
	
	public:

	    /* Type definitions: */

	    typedef typename base_type::size_type size_type;
	    typedef value_t value_type;
	    typedef typename base_type::object_type object_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    
	    /// Area of the D-dimensional unit sphere.

	    static const value_t unit_sphere_area;

	    /// Radius address.

	    const value_type* radius;

	    /// Default constructor.

	    uniform_sphere(const value_type* radius_=NULL):radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::pow(std::abs(*radius),(int)D));
		this->weight()=vol;
	    }

	    /// Constructor with point address argument.

	    uniform_sphere(object_type* instance_,const value_type* radius_=NULL):base_type(instance_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::pow(std::abs(*radius),(int)D));
		this->weight()=vol;
	    }

	    /// Constructor with point address and integrand address arguments.

	    uniform_sphere(object_type* instance_,const value_type* integrand_,const value_type* radius_=NULL):base_type(instance_,integrand_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::pow(std::abs(*radius),(int)D));
		this->weight()=vol;
	    }

	    /// Clone method implementation

	    uniform_sphere<value_t,D,rng_t>* clone() const
	    {
		return new uniform_sphere<value_t,D,rng_t>(*this);
	    }

	    /// Generation method.

	    bool generate()
	    {
		value_type w=(value_type)0;
		for(size_type i=0;i<D+1;++i)
		{
		    Gauss_gen->generate();
		    (this->object())[i]=x;
		    w+=(x*x);
		}
		w=(radius==NULL)?((value_type)1/std::sqrt(w)):(std::abs(*radius)/std::sqrt(w));
		(this->object())*=w;
		this->weight()=vol;
		return true;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		this->weight()=vol;
		return true;
	    }

	    /// Parameter refresher function.

	    void refresh_params()
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::pow(std::abs(*radius),(int)D));
	    }
	    
	private:

	    value_type vol;
	    static value_t x;
	    static normal_generator<value_t,rng_t>* Gauss_gen;
    };
    template<class value_t,std::size_t D,class rng_t>const value_t uniform_sphere<value_t,D,rng_t>::unit_sphere_area=(D%2==0)?((value_t)2*std::pow(std::acos(-(value_t)1),int(D/2))/double_factorial(D-1)):(std::pow(std::acos(-(value_t)1),int((D+1)/2))/double_factorial(D-1));
    template<class value_t,std::size_t D,class rng_t>value_t uniform_sphere<value_t,D,rng_t>::x(0);
    template<class value_t,std::size_t D,class rng_t>normal_generator<value_t,rng_t>* uniform_sphere<value_t,D,rng_t>::Gauss_gen=new normal_generator<value_t,rng_t>(&uniform_sphere<value_t,D,rng_t>::x);
    
    template<class value_t,class rng_t>class uniform_sphere<value_t,1,rng_t>: public object_allocator< vector<value_t,2> >, public MC_generator<value_t>
    {
	typedef object_allocator< vector<value_t,2> > base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::object_type object_type;
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    
	    /// Length of the unit circle.

	    static const value_t unit_sphere_area;

	    /// Radius address.

	    const value_type* radius;

	    /// Default constructor.

	    uniform_sphere(const value_type* radius_=NULL):radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::abs(*radius));
		this->weight()=vol;
	    }

	    /// Constructor with point address argument.

	    uniform_sphere(object_type* instance_,const value_type* radius_=NULL):base_type(instance_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::abs(*radius));
		this->weight()=vol;
	    }

	    /// Constructor with point address and integrand address arguments.

	    uniform_sphere(object_type* instance_,const value_type* integrand_,const value_type* radius_=NULL):base_type(instance_,integrand_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::abs(*radius));
		this->weight()=vol;
	    }

	    /// Clone method implementation

	    uniform_sphere<value_t,1,rng_t>* clone() const
	    {
		return new uniform_sphere<value_t,1,rng_t>(*this);
	    }

	    /// Generation method.

	    bool generate()
	    {
		value_type w,x0,x1;
		do
		{
		    x0=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    x1=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    w=x0*x0+x1*x1;
		}
		while(w>=(value_type)1);
		(this->object())[0]=(radius==NULL)?((x0*x0-x1*x1)/w):((*radius)*(x0*x0-x1*x1)/w);
		(this->object())[1]=(radius==NULL)?((value_type)2*x0*x1/w):((value_type)2*(*radius)*x0*x1/w);
		this->weight()=vol;
		return true;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		this->weight()=vol;
		return true;
	    }

	    /// Parameter refresher function.

	    void refresh_params()
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*std::abs(*radius));
	    }
	    
	private:

	    value_type vol;
    };
    template<class value_t,class rng_t>const value_t uniform_sphere<value_t,1,rng_t>::unit_sphere_area=(value_t)2*std::acos(-(value_t)1);
    
    template<class value_t,class rng_t>class uniform_sphere<value_t,2,rng_t>: public object_allocator< vector<value_t,3> >, public MC_generator<value_t>
    {
	typedef object_allocator< vector<value_t,3> > base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::object_type object_type;
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    
	    /// Area of the 2-dimensional unit sphere.

	    static const value_t unit_sphere_area;

	    /// Radius address.

	    const value_type* radius;

	    /// Default constructor.

	    uniform_sphere(const value_type* radius_=NULL):radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius));
		this->weight()=vol;
	    }

	    /// Constructor with point address argument.

	    uniform_sphere(object_type* instance_,const value_type* radius_=NULL):base_type(instance_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius));
		this->weight()=vol;
	    }

	    /// Constructor with point address and integrand address arguments.

	    uniform_sphere(object_type* instance_,const value_type* integrand_,const value_type* radius_=NULL):base_type(instance_,integrand_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius));
		this->weight()=vol;
	    }

	    /// Clone method implementation

	    uniform_sphere<value_t,2,rng_t>* clone() const
	    {
		return new uniform_sphere<value_t,2,rng_t>(*this);
	    }

	    /// Generation method.

	    bool generate()
	    {
		value_type w,x0,x1;
		do
		{
		    x0=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    x1=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    w=x0*x0+x1*x1;
		}
		while(w>=(value_type)1);
		value_type y=(radius==NULL)?((value_type)2*std::sqrt((value_type)1-w)):((value_type)2*(*radius)*std::sqrt((value_type)1-w));
		(this->object())[0]=y*x0;
		(this->object())[1]=y*x1;
		(this->object())[2]=(radius==NULL)?((value_type)1-(value_type)2*w):(*radius-(value_t)2*(*radius)*w);
		this->weight()=vol;
		return true;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		this->weight()=vol;
		return true;
	    }

	    /// Parameter refresher function.

	    void refresh_params()
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius));
	    }
	    
	private:

	    value_type vol;
    };
    template<class value_t,class rng_t>const value_t uniform_sphere<value_t,2,rng_t>::unit_sphere_area=(value_t)4*std::acos(-(value_t)1);
    
    template<class value_t,class rng_t>class uniform_sphere<value_t,3,rng_t>: public object_allocator< vector<value_t,4> >, public MC_generator<value_t>
    {
	typedef object_allocator< vector<value_t,4> > base_type;

	public:

	    /* Type definitions: */

	    typedef typename base_type::size_type size_type;
	    typedef typename base_type::object_type object_type;
	    typedef value_t value_type;
	    typedef rng_t rn_engine;
	    typedef random_number_stream<value_t,rng_t> rn_stream;
	    
	    /// Area of the 3-dimensional unit hypersphere.

	    static const value_t unit_sphere_area;

	    /// Radius address.

	    const value_type* radius;

	    /// Default constructor.

	    uniform_sphere(const value_type* radius_=NULL):radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius))*std::abs(*radius);
		this->weight()=vol;
	    }

	    /// Constructor with point address argument.

	    uniform_sphere(object_type* instance_,const value_type* radius_=NULL):base_type(instance_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius))*std::abs(*radius);
		this->weight()=vol;
	    }

	    /// Constructor with point address and integrand address arguments.

	    uniform_sphere(object_type* instance_,const value_type* integrand_,const value_type* radius_=NULL):base_type(instance_,integrand_),radius(radius_)
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius)*std::abs(*radius));
		this->weight()=vol;
	    }

	    /// Clone method implementation

	    uniform_sphere<value_t,3,rng_t>* clone() const
	    {
		return new uniform_sphere<value_t,3,rng_t>(*this);
	    }

	    /// Generation method.

	    bool generate()
	    {
		value_type w0,w1,x0,x1,x2,x3;
		do
		{
		    x0=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    x1=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    x2=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    x3=rn_stream::throw_number(-(value_type)1,(value_type)1);
		    w0=x0*x0+x1*x1;
		    w1=x2*x2+x3*x3;
		}
		while(w0>=(value_t)1 or w1>=(value_t)1);
		(this->object())[0]=(radius==NULL)?x0:((*radius)*x0);
		(this->object())[1]=(radius==NULL)?x1:((*radius)*x1);
		value_t z=std::sqrt(((value_t)1-w0)/w1);
		(this->object())[2]=(radius==NULL)?(z*x2):(z*(*radius)*x2);
		(this->object())[3]=(radius==NULL)?(z*x3):(z*(*radius)*x3);
		this->weight()=vol;
		return true;
	    }

	    /// Weight evaluation method.

	    bool evaluate_weight()
	    {
		this->weight()=vol;
		return true;
	    }

	    /// Parameter refresher function.

	    void refresh_params()
	    {
		vol=(radius==NULL)?(unit_sphere_area):(unit_sphere_area*(*radius)*(*radius))*std::abs(*radius);
	    }
	    
	private:

	    value_type vol;
    };
    template<class value_t,class rng_t>const value_t uniform_sphere<value_t,3,rng_t>::unit_sphere_area=(value_t)2*std::pow(std::acos(-(value_t)1),2);
}

#endif /*CAMGEN_UNI_SPHERE_H_*/

