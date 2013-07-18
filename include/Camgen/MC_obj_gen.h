//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file MC_obj_gen.h
    \brief Base class template for Monte Carlo generator classes with object address member.
 */

#ifndef CAMGEN_MC_OBJ_GEN_H_
#define CAMGEN_MC_OBJ_GEN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Adds extra functionality to the Monte Carlo generator classes by including an *
 * object address data member to be filled by the derived class' generate method *
 * implementation.                                                               *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */                                                                              
#include <bitset>
#include <Camgen/vector.h>
#include <Camgen/MC_gen.h>

namespace Camgen
{
    /// Object instance containing MC generator class. The first template parameter
    /// denotes the numerical type of the weights, the second parameter the object instance
    /// type, the third argument the number of instances to be handled and the fourth
    /// parameter a flag denoting whether the object type contains a default constructor.

    template<class value_t,class object_t,std::size_t N=1,bool q=true>class MC_object_generator;

    /// Specialization for objects with default constructor.

    template<class value_t,class object_t,std::size_t N>class MC_object_generator<value_t,object_t,N,true>: public MC_generator<value_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef std::size_t size_type;
	    typedef object_t object_type;
	    typedef typename MC_generator<value_t>::integral_type integral_type;

	    /// Default constructor, allocating the instances with default constructors.

	    MC_object_generator()
	    {
		for(size_type i=0;i<N;++i)
		{
		    instances[i]=new object_type;
		    alloc_obj[i]=true;
		}
	    }

	    /// Constructor with input instances.

	    MC_object_generator(const vector<object_type*,N>& instances_):instances(instances_)
	    {
		alloc_obj.reset();
	    }

	    /// Copy constructor.

	    MC_object_generator(const MC_object_generator<value_t,object_t,N,true>& other):MC_generator<value_t>(other),alloc_obj(other.alloc_obj)
	    {
		for(size_type i=0;i<N;++i)
		{
		    if(alloc_obj[i])
		    {
			instances[i]=new object_type(*(other.instances[i]));
		    }
		}
	    }

	    /// Destructor.

	    virtual ~MC_object_generator()
	    {
		for(size_type i=0;i<N;++i)
		{
		    if(alloc_obj[i])
		    {
			delete instances[i];
		    }
		}
	    }

	    /// Returns the reference to the i-th instance.
	    
	    object_type& object(size_type i=0)
	    {
		return *(instances[i]);
	    }

	    /// Returns the const reference to the i-th instance.
	    
	    const object_type& object(size_type i=0) const
	    {
		return *(instances[i]);
	    }

	    /// Returns the address of the i-th instance.
	    
	    object_type* get_object(size_type i=0)
	    {
		return instances[i];
	    }

	    /// Returns the const pointer to the i-th instance.
	    
	    const object_type* get_object(size_type i=0) const
	    {
		return instances[i];
	    }

	    /// Sets the address to the i-th instance.

	    void set_object(object_type* instance_,size_type i=0)
	    {
		if(alloc_obj[i])
		{
		    delete instances[i];
		    alloc_obj.reset(i);
		}
		instances[i]=instance_;
	    }

	    /// Calls the generation method and returns the i-th object.

	    const object_type& operator ()(size_type i)
	    {
		this->generate();
		return *(instances[i]);
	    }

	private:

	    vector<object_type*,N> instances;
	    std::bitset<N> alloc_obj;
    };

    /// Specialization for objects without default constructor.

    template<class value_t,class object_t,std::size_t N>class MC_object_generator<value_t,object_t,N,false>: public MC_generator<value_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef std::size_t size_type;
	    typedef object_t object_type;
	    typedef typename MC_generator<value_t>::integral_type integral_type;

	    /// Constructor with input instances.

	    MC_object_generator(const vector<object_type*,N>& instances_):instances(instances_){}

	    /// Copy constructor.

	    MC_object_generator(const MC_object_generator<value_t,object_t,N,false>& other):MC_generator<value_t>(other),instances(other.instances){}

	    /// Destructor.

	    virtual ~MC_object_generator(){}

	    /// Returns the reference to the i-th instance.
	    
	    object_type& object(size_type i=0)
	    {
		return *(instances[i]);
	    }

	    /// Returns the const reference to the i-th instance.
	    
	    const object_type& object(size_type i=0) const
	    {
		return *(instances[i]);
	    }

	    /// Returns the address of the i-th instance.
	    
	    object_type* get_object(size_type i=0)
	    {
		return instances[i];
	    }

	    /// Returns the const pointer to the i-th instance.
	    
	    const object_type* get_object(size_type i=0) const
	    {
		return instances[i];
	    }

	    /// Sets the address to the i-th instance.

	    void set_object(object_type* instance_,size_type i=0)
	    {
		instances[i]=instance_;
	    }

	    /// Calls the generation method and returns the i-th object.

	    const object_type& operator ()(size_type i)
	    {
		this->generate();
		return *(instances[i]);
	    }

	private:

	    vector<object_type*,N> instances;
    };

    /// Specialization for objects with default constructor.

    template<class value_t,class object_t>class MC_object_generator<value_t,object_t,1,true>: public MC_generator<value_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef std::size_t size_type;
	    typedef object_t object_type;
	    typedef typename MC_generator<value_t>::integral_type integral_type;

	    /// Default constructor, allocating the instances with default constructors.

	    MC_object_generator():instance(new object_type),alloc_obj(true){}

	    /// Constructor with input instance.

	    MC_object_generator(object_type* instance_):instance(instance_),alloc_obj(false){}

	    /// Copy constructor.

	    MC_object_generator(const MC_object_generator<value_t,object_t,1,true>& other):MC_generator<value_t>(other),alloc_obj(other.alloc_obj)
	    {
		if(alloc_obj)
		{
		    instance=new object_type(*(other.instance));
		}
		else
		{
		    instance=other.instance;
		}
	    }

	    /// Destructor.

	    virtual ~MC_object_generator()
	    {
		if(alloc_obj)
		{
		    delete instance;
		}
	    }

	    /// Returns the reference to the i-th instance.
	    
	    object_type& object()
	    {
		return *instance;
	    }

	    /// Returns the const reference to the i-th instance.
	    
	    const object_type& object() const
	    {
		return *instance;
	    }

	    /// Returns the address of the i-th instance.
	    
	    object_type* get_object()
	    {
		return instance;
	    }

	    /// Returns the const pointer to the i-th instance.
	    
	    const object_type* get_object() const
	    {
		return instance;
	    }

	    /// Sets the address to the i-th instance.

	    void set_object(object_type* instance_)
	    {
		if(alloc_obj)
		{
		    delete instance;
		    alloc_obj=false;
		}
		instance=instance_;
	    }

	    /// Calls the generation method and returns the object.

	    const object_type& operator ()()
	    {
		this->generate();
		return *instance;
	    }

	private:

	    object_type* instance;
	    bool alloc_obj;
    };

    /// Specialization for objects without default constructor.

    template<class value_t,class object_t>class MC_object_generator<value_t,object_t,1,false>: public MC_generator<value_t>
    {
	public:

	    /* Type definitions: */

	    typedef value_t value_type;
	    typedef std::size_t size_type;
	    typedef object_t object_type;
	    typedef typename MC_generator<value_t>::integral_type integral_type;

	    /// Constructor with input instance.

	    MC_object_generator(object_type* instance_):instance(instance_){}

	    /// Copy constructor.

	    MC_object_generator(const MC_object_generator<value_t,object_t,1,false>& other):MC_generator<value_t>(other),instance(other.instance){}

	    /// Destructor.

	    virtual ~MC_object_generator(){}

	    /// Returns the reference to the i-th instance.
	    
	    object_type& object()
	    {
		return *instance;
	    }

	    /// Returns the const reference to the i-th instance.
	    
	    const object_type& object() const
	    {
		return *instance;
	    }

	    /// Returns the address of the i-th instance.
	    
	    object_type* get_object()
	    {
		return instance;
	    }

	    /// Returns the const pointer to the i-th instance.
	    
	    const object_type* get_object() const
	    {
		return instance;
	    }

	    /// Sets the address to the i-th instance.

	    void set_object(object_type* instance_)
	    {
		instance=instance_;
	    }

	    /// Calls the generation method and returns the object.

	    const object_type& operator ()()
	    {
		this->generate();
		return *instance;
	    }

	private:

	    object_type* instance;
    };
}

#endif /*CAMGEN_MC_OBJ_GEN_H_*/

