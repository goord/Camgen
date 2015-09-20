//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file obj_alloc.h
    \brief Base class template for allocating member objects.
 */

#ifndef CAMGEN_OBJ_ALLOC_H_
#define CAMGEN_OBJ_ALLOC_H_

#include <vector>

namespace Camgen
{
    /// Class that automatically allocates a series of objects of type object_t, assuming these contain default
    /// constructors.

    template<class object_t>class object_allocator
    {
	public:

	    /* Type definitions: */

	    typedef object_t object_type;
	    typedef std::vector<object_type*> instance_vector_type;
	    typedef typename instance_vector_type::size_type size_type;

	    /// Constructor allocating n instances of object_type.

	    object_allocator(size_type n=1):alloc_obj(n,true)
	    {
		for(size_type i=0;i<n;++i)
		{
		    instances.push_back(new object_type);
		}
	    }

	    /// Non-allocating constructor, copying the instances from the argument vector.

	    object_allocator(const std::vector<object_type*>& instances_):instances(instances_),alloc_obj(instances_.size(),false){}

	    /// Non-allocating constructor with one object instance.

	    object_allocator(object_type* instance_):instances(1,instance_),alloc_obj(1,false){}

	    /// Non-allocating constructor with vector iterator arguments.
	    
	    template<class input_iterator> object_allocator(input_iterator start,input_iterator end):instances(start,end),alloc_obj(instances.size(),false){}

	    /// Copy constructor.

	    object_allocator(const object_allocator<object_t>& other):alloc_obj(other.alloc_obj)
	    {
		for(size_type i=0;i<alloc_obj.size();++i)
		{
		    instances.push_back(alloc_obj[i]?new object_type(*(other.instances[i])):other.instances[i]);
		}
	    }

	    /// Destructor.

	    virtual ~object_allocator()
	    {
		for(size_type i=0;i<size();++i)
		{
		    if(alloc_obj[i])
		    {
			delete instances[i];
		    }
		}
	    }

	    /* Public readout methods: */
	    /*-------------------------*/

	    /// Returns the number on instances: */
	    
	    size_type size() const
	    {
		return instances.size();
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

	    /* Public modifier methods: */
	    /*--------------------------*/

	    /// Sets the address to the i-th instance.

	    void set_object(object_type* instance_,size_type i=0)
	    {
		if(alloc_obj[i])
		{
		    delete instances[i];
		    alloc_obj[i]=false;
		}
		instances[i]=instance_;
	    }

	private:

	    std::vector<object_type*> instances;
	    std::vector<bool> alloc_obj;
    };
}

#endif /*CAMGEN_OBJ_ALLOC_H_*/

