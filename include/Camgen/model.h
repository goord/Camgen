//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file model.h
    \brief model base class template interface and implementation header.
 */

#ifndef CAMGEN_MODEL_H_
#define CAMGEN_MODEL_H_

#include <ostream>
#include <fstream>
#include <string>
#include <cxxabi.h>
#include <Camgen/unpack_vert.h>
#include <Camgen/license_print.h>
#include <Camgen/model_wrapper.h>
#include <Camgen/scalar_particle.h>
#include <Camgen/fermion.h>
#include <Camgen/vector_particle.h>
#include <Camgen/comp_vert.h>
#include <Camgen/vertex.h>
#include <Camgen/adj_rep.h>
#include <Camgen/fundam_rep.h>
#include <Camgen/f.h>
#include <Camgen/asymtvv.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Model base class definition. The model base class takes the implementation    *
 * type as template parameter, and assumes that the constructor of this          *
 * implementation class consists of a series of statements adding particles and  *
 * vertices (in that order). The base class has a singleton structure, and as    *
 * such the implementation constructor shall be called only once. Furthermore it *
 * outputs all actions involving adding or erasing of particles to a separate    *
 * log file called "model.log", useful for debugging a model.                    *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    /// model base class template.
    /** All user-model classes should be derived from model<derived-model-type>.
     * This class has no data members or non-static member functions. It serves
     * only as an interface between the derived model_t implementation and
     * the hidden model_wrapper<model_t> class templte, which contains the vertex and particle content. */

    template<class model_t>class model
    {
	public:

	    /// Returns the single instance of model_t.
	    /** If the instance is not yet allocated, uses the default
	     * constructor of model_t to create and return it. */

	    static model_t* get_instance()
	    {
		license_print::initialise();
		if(instance==NULL)
		{
		    instance=new model_t;
		}
		return instance;
	    }

	    /// Creates the single instance of model_t.
	    /** Same functionality as the get_instance(), but no return value.
	     * */ 

	    static void initialise()
	    {
		license_print::initialise();
		if(instance==NULL)
		{
		    instance=new model_t;
		}
	    }

	    /// Requests whether the singleton instance was created.

	    static bool initialised()
	    {
		return (obj_counter > 0);
	    }

	    /// Returns the log file stream.

	    static std::ostream& logfile()
	    {
		return fs;
	    }

	    /// Inserts spin-s auxiliary particle with name str.
	    /// Auxiliary particles in Camgen do not propagate and have no wave
	    /// functions. They can therefore not be selected as external
	    /// particles.

	    static void add_auxiliary_particle(const std::string& str,const spin& s)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::create_auxiliary(str,s));
		}
	    }

	    /// Inserts a pair of conjugate spin-s auxiliary particles with names (str1,str2).

	    static void add_auxiliary_particles(const std::string& str1,const std::string& str2,const spin& s)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::create_auxiliary(str1,str2,s));
		}
	    }

	    /// Inserts spin-s auxiliary multiplet with name str in representation rep_t.
	    /// The rep_t should be either adjoint_rep<group_t> or
	    /// fundamental_rep<group_t> for some group type group_t, or a
	    /// composition thereof using the compose_reps class template. Currently
	    /// only the SU<N> classes may serve as a group type.

	    template<class rep_t>static void add_auxiliary_particle(const std::string& str,const spin& s)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create_auxiliary<rep_t>(str,s));
		}
	    }

	    /// Inserts a pair of conjugate spin-s auxiliary multiplets with names (str1,str2) in representation rep_t.

	    template<class rep_t>static void add_auxiliary_particles(const std::string& str1,const std::string& str2,const spin& s)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create_auxiliary<rep_t>(str1,str2,s));
		}
	    }

	    /// Inserts a massless scalar particle with name str and PDG-id n.

	    static void add_scalar(const std::string& str,int n=0)
	    {
		if(obj_counter<2)
		{	
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create< scalar_particle<model_t> >(str,NULL,NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate massless scalar particles with names (str1,str2) and PDG-ids (n,-n).
	    
	    static void add_scalars(const std::string& str1,const std::string str2,int n=0)
	    {
		if(obj_counter<2)
		{	
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create< scalar_particle<model_t> >(str1,str2,NULL,NULL,n));
		}
	    }

	    /// Inserts a scalar particle with name str, mass *m and PDG-id n.

	    static void add_scalar(const std::string& str,void* m,int n=0)
	    {
		if(obj_counter<2)
		{	
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create< scalar_particle<model_t> >(str,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate scalar particles with names (str1,str2), mass *m and PDG-ids (n,-n).
	    
	    static void add_scalars(const std::string& str1,const std::string& str2,void* m,int n=0)
	    {
		if(obj_counter<2)
		{	
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create< scalar_particle<model_t> >(str1,str2,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a scalar particle with name str, mass *m, width *w and PDG-id n.

	    static void add_scalar(const std::string& str,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{	
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create< scalar_particle<model_t> >(str,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a pair of conjugate scalar particles with names (str1,str2), mass *m, width *w and PDG-ids (n,-n).

	    static void add_scalars(const std::string& str1,const std::string& str2,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{	
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create< scalar_particle<model_t> >(str1,str2,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a massless scalar multiplet with name str and PDG-id n in representation rep_t.

	    template<class rep_t>static void add_scalar(const std::string& str,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<scalar_particle<model_t>,rep_t>(str,NULL,NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate massless scalar multiplets with names (str1,str2) and PDG-ids (n,-n) in representation rep_t.
	    
	    template<class rep_t>static void add_scalars(const std::string& str1,const std::string& str2,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<scalar_particle<model_t>,rep_t>(str1,str2,NULL,NULL,n));
		}
	    }

	    /// Inserts a scalar multiplet with name str, mass *m and PDG-id n in representation rep_t.

	    template<class rep_t>static void add_scalar(const std::string& str,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<scalar_particle<model_t>,rep_t>(str,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate scalar multiplets with names (str1,str2), mass *m and PDG-ids (n,-n) in representation rep_t.
	    
	    template<class rep_t>static void add_scalars(const std::string& str1,const std::string& str2,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<scalar_particle<model_t>,rep_t>(str1,str2,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a scalar multiplet with name str, mass *m, width *w and PDG-id n in representation rep_t.

	    template<class rep_t>static void add_scalar(const std::string& str,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<scalar_particle<model_t>,rep_t>(str,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a pair of conjugate scalar multiplets with names (str1,str2), mass *m, width *w and PDG-ids (n,-n) in representation rep_t.
	    
	    template<class rep_t>static void add_scalars(const std::string& str1,const std::string& str2,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<scalar_particle<model_t>,rep_t>(str1,str2,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a massless Majorana fermion with name str and PDG-id n.

	    static void add_fermion(const std::string& str,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create< fermion<model_t,typename model_t::spacetime_type,model_t::dimension> >(str,NULL,NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate Dirac fermions with names (str1,str2) and PDG-ids (n,-n).
	    /// The first argument denotes the particle and the second one the
	    /// anti-particle, e.g. add_fermions("e-","e+",11) in QED.
	    
	    static void add_fermions(const std::string& str1,const std::string& str2,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create< fermion<model_t,typename model_t::spacetime_type,model_t::dimension> >(str1,str2,NULL,NULL,n));
		}
	    }

	    /// Inserts a Majorana fermion with name str, mass *m and PDG-id n.

	    static void add_fermion(const std::string& str,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create< fermion<model_t,typename model_t::spacetime_type,model_t::dimension> >(str,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate Dirac fermions with names (str1,str2), mass *m and PDG-ids (n,-n).

	    static void add_fermions(const std::string& str1,const std::string& str2,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create< fermion<model_t,typename model_t::spacetime_type,model_t::dimension> >(str1,str2,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a Majorana fermion with name str, mass *m, width *w and PDG-id n.

	    static void add_fermion(const std::string& str,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create< fermion<model_t,typename model_t::spacetime_type,model_t::dimension> >(str,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a pair of conjugate Dirac fermions with names (str1,str2), mass *m, width *w and PDG-ids (n,-n).

	    static void add_fermions(const std::string& str1,const std::string& str2,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create< fermion<model_t,typename model_t::spacetime_type,model_t::dimension> >(str1,str2,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a massless Majorana fermion multiplet with name str and PDG-id n in representation rep_t.

	    template<class rep_t>static void add_fermion(const std::string& str,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<fermion<model_t,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str,NULL,NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate Dirac fermion multiplets with names (str1,str2) and PDG-ids (n,-n) in representation rep_t.
	    
	    template<class rep_t>static void add_fermions(const std::string& str1,const std::string& str2,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<fermion<model_t,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str1,str2,NULL,NULL,n));
		}
	    }

	    /// Inserts a Majorana fermion with name str, mass *m and PDG-id n in representation rep_t.

	    template<class rep_t>static void add_fermion(const std::string& str,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<fermion<model_t,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate Dirac fermions with names (str1,str2), mass *m and PDG-ids (n,-n) in representation rep_t.

	    template<class rep_t>static void add_fermions(const std::string& str1,const std::string& str2,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<fermion<model_t,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str1,str2,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a Majorana fermion with name str, mass *m, width *w and PDG-id n in representation rep_t.

	    template<class rep_t>static void add_fermion(const std::string& str,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<fermion<model_t,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a pair of conjugate Dirac fermions with names (str1,str2), mass *m, width *w and PDG-ids (n,-n) in representation rep_t.

	    template<class rep_t>static void add_fermions(const std::string& str1,const std::string& str2,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<fermion<model_t,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str1,str2,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a massless vector particle in gauge g with name str and PDG-id n.

	    template<template<class mod>class g>static void add_vector(const std::string& str,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension> >(str,NULL,NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate massless vector particles with names (str1,str2) and PDG-ids (n,-n).

	    template<template<class mod>class g>static void add_vectors(const std::string& str1,const std::string& str2,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension> >(str1,str2,NULL,NULL,n));
		}
	    }

	    /// Inserts a vector particle in gauge g with name str, mass *m and PDG-id n.

	    template<template<class mod>class g>static void add_vector(const std::string& str,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension> >(str,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate vector particles with names (str1,str2), mass *m and PDG-ids (n,-n).

	    template<template<class mod>class g>static void add_vectors(const std::string& str1,const std::string& str2,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension> >(str1,str2,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a vector particle in gauge g with name str, mass *m, width *w and PDG-id n.

	    template<template<class mod>class g>static void add_vector(const std::string& str,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension> >(str,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a pair of conjugate vector particles with names (str1,str2), mass *m, width *w and PDG-ids (n,-n).

	    template<template<class mod>class g>static void add_vectors(const std::string& str1,const std::string& str2,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension> >(str1,str2,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a massless vector multiplet in gauge g with name str and PDG-id n in representation rep_t.

	    template<class rep_t,template<class mod>class g>static void add_vector(const std::string& str,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str,NULL,NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate massless vector multiplets with names (str1,str2) and PDG-ids (n,-n) in representation rep_t.

	    template<class rep_t,template<class mod>class g>static void add_vectors(const std::string& str1,const std::string& str2,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str1,str2,NULL,NULL,n));
		}
	    }

	    /// Inserts a vector multiplet in gauge g with name str, mass *m and PDG-id n in representation rep_t.

	    template<class rep_t,template<class mod>class g>static void add_vector(const std::string& str,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a pair of conjugate vector multiplets with names (str1,str2), mass *m and PDG-ids (n,-n) in representation rep_t.

	    template<class rep_t,template<class mod>class g>static void add_vectors(const std::string& str1,const std::string& str2,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str1,str2,static_cast<const typename model_t::value_type*>(m),NULL,n));
		}
	    }

	    /// Inserts a vector multiplet in gauge g with name str, mass *m, width *w and PDG-id n in representation rep_t.

	    template<class rep_t,template<class mod>class g>static void add_vector(const std::string& str,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Inserts a pair of conjugate vector multiplets with names (str1,str2), mass *m, width *w and PDG-ids (n,-n) in representation rep_t.

	    template<class rep_t,template<class mod>class g>static void add_vectors(const std::string& str1,const std::string& str2,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_particle_pair(particle<model_t>::template create<vector_particle<model_t,g<model_t>,typename model_t::spacetime_type,model_t::dimension>,rep_t>(str1,str2,static_cast<const typename model_t::value_type*>(m),static_cast<const typename model_t::value_type*>(w),n));
		}
	    }

	    /// Insert a photon in gauge g (name "gamma", PDG-id 22).
	    /// Equivalent to: add_vector<g>("gamma",22).

	    template<template<class mod>class g>static void add_photon()
	    {
		add_vector<g>("gamma",22);
	    }

	    /// Insert a gluon of group group_t in gauge g (name "g", PDG-id 21).
	    /// Equivalent to: add_vector<adjoint_rep<group_t>,g>("g",21).

	    template<template<class mod>class g,class group_t>static void add_gluon()
	    {
		if(obj_counter<2)
		{
		    add_vector<adjoint_rep<group_t>,g>("g",21);
		}
	    }

	    /// Inserts massless quarks with names (str1,str2) and PDG-ids (n,-n) in the fundamental representation of group_t.
	    /// Equivalent to: add_fermions< fundamental_rep<group_t>(str1,str2,n).

	    template<class group_t>static void add_quarks(const std::string& str1,const std::string& str2,int n=0)
	    {
		if(obj_counter<2)
		{
		    add_fermions< fundamental_rep<group_t> >(str1,str2,n);
		}
	    }

	    /// Inserts quarks with names (str1,str2), mass *m and PDG-ids (n,-n) in the fundamental representation of group_t.
	    /// Equivalent to: add_fermions< fundamental_rep<group_t>(str1,str2,m,n).

	    template<class group_t>static void add_quarks(const std::string& str1,const std::string str2,void* m,int n=0)
	    {
		if(obj_counter<2)
		{
		    add_fermions< fundamental_rep<group_t> >(str1,str2,m,n);
		}
	    }

	    /// Inserts quarks with names (str1,str2), mass *m, width *w and PDG-ids (n,-n) in the fundamental representation of group_t.
	    /// Equivalent to: add_fermions< fundamental_rep<group_t>(str1,str2,m,n).

	    template<class group_t>static void add_quarks(const std::string& str1,const std::string str2,void* m,void* w,int n=0)
	    {
		if(obj_counter<2)
		{
		    add_fermions< fundamental_rep<group_t> >(str1,str2,m,w,n);
		}
	    }

	    /// Inserts a 3-vertex connecting particles phi1,phi2,phi3 by Feynman rule Feynrule_t with couplings *c1(,*c2,*c3).
	    /// Note that the couplings in vertex insertions should be of the type std::complex<value_type>*

	    template<template<class mod>class Feynrule_t>static void add_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3,void* c1,void* c2=NULL,void* c3=NULL)
	    {
		const std::complex<typename model_t::value_type>* d2=(c2==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c2);
		const std::complex<typename model_t::value_type>* d3=(c3==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c3);

		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_vertex(vertex<model_t>::template create<Feynrule_t<model_t> >(model_wrapper<model_t>::get_particle(phi1),model_wrapper<model_t>::get_particle(phi2),model_wrapper<model_t>::get_particle(phi3),static_cast<const std::complex<typename model_t::value_type>*>(c1),d2,d3));
		}
	    }

	    /// Inserts a 3-vertex connecting particles phi1,phi2,phi3 by Feynman rule Feynrule_t with colour structure coltens_t and couplings *c1(,*c2,*c3).

	    template<class coltens_t,template<class mod>class Feynrule_t>static void add_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3,void* c1,void* c2=NULL,void* c3=NULL)
	    {
		const std::complex<typename model_t::value_type>* d2=(c2==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c2);
		const std::complex<typename model_t::value_type>* d3=(c3==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c3);
		
		if(obj_counter<2)
		{
		    typedef typename unpack_vertices< coltens_t,Feynrule_t<model_t> >::type inner_type;
		    typedef typename get_colour_treatment<model_t,model_t::coloured>::type::template make_Feynrule<inner_type>::type vert_type;

		    model_wrapper<model_t>::insert_vertex(vertex<model_t>::template create<vert_type>(model_wrapper<model_t>::get_particle(phi1),model_wrapper<model_t>::get_particle(phi2),model_wrapper<model_t>::get_particle(phi3),static_cast<const std::complex<typename model_t::value_type>*>(c1),d2,d3));
		}
	    }

	    /// Inserts a 4-vertex connecting particles phi1,phi2,phi3,phi4 by Feynman rule Feynrule_t with couplings *c1(,*c2,*c3).

	    template<template<class mod>class Feynrule_t>static void add_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3,const std::string& phi4,void* c1,void* c2=NULL,void* c3=NULL)
	    {
		const std::complex<typename model_t::value_type>* d2=(c2==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c2);
		const std::complex<typename model_t::value_type>* d3=(c3==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c3);

		if(obj_counter<2)
		{
		    model_wrapper<model_t>::insert_vertex(vertex<model_t>::template create<Feynrule_t<model_t> >(model_wrapper<model_t>::get_particle(phi1),model_wrapper<model_t>::get_particle(phi2),model_wrapper<model_t>::get_particle(phi3),model_wrapper<model_t>::get_particle(phi4),static_cast<const std::complex<typename model_t::value_type>*>(c1),d2,d3));
		}
	    }

	    /// Inserts a 4-vertex connecting particles phi1,phi2,phi3,phi4 by Feynman rule Feynrule_t with colour structure coltens_t and couplings *c1(,*c2,*c3).

	    template<class coltens_t,template<class mod>class Feynrule_t>static void add_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3,const std::string& phi4,void* c1,void* c2=NULL,void* c3=NULL)
	    {
		const std::complex<typename model_t::value_type>* d2=(c2==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c2);
		const std::complex<typename model_t::value_type>* d3=(c3==NULL)?NULL:static_cast<const std::complex<typename model_t::value_type>*>(c3);
		
		if(obj_counter<2)
		{
		    typedef typename unpack_vertices< coltens_t,Feynrule_t<model_t> >::type inner_type;
		    typedef typename get_colour_treatment<model_t,model_t::coloured>::type::template make_Feynrule<inner_type>::type vert_type;

		    model_wrapper<model_t>::insert_vertex(vertex<model_t>::template create<vert_type>(model_wrapper<model_t>::get_particle(phi1),model_wrapper<model_t>::get_particle(phi2),model_wrapper<model_t>::get_particle(phi3),model_wrapper<model_t>::get_particle(phi4),static_cast<const std::complex<typename model_t::value_type>*>(c1),d2,d3));
		}
	    }

	    /// Insert the fast 4-gluon vertex.
	    /// Inserts an antisymmetric tensor field in the
	    /// adjoint representation and a vertex of coupling it to 2 phi-particles
	    /// an antisymmetric tensor coupling composed with f< SU<N> >. This will optimize gluon
	    /// process performance.

	    template<class group_t>static void add_fast_4g_vertex(const std::string& phi,void* c0)
	    {
		add_auxiliary_particle< adjoint_rep<group_t> >("Hqcd",2);
		add_vertex<colour_tensor::f<group_t,0,1,2>,asymtvv>("Hqcd",phi,phi,c0,NULL,NULL);
	    }

	    /// Inserts a particle family with name str and content cont to the model wrapper.

	    static void construct_family(const std::string& str,const std::string& cont)
	    {
		if(obj_counter<2)
		{
		    model_wrapper<model_t>::construct_family(str,cont);
		}
	    }

	    /// Sends compile-time data, particle content and vertex list to os.

	    static std::ostream& print(std::ostream& os)
	    {
		os<<"dimension:               "<<model_t::dimension<<std::endl;
		os<<std::endl;
		typename model_t::value_type x;
		os<<"value type:              "<<abi::__cxa_demangle(typeid(x).name(),0,0,NULL)<<std::endl;
		os<<std::endl;
		typename get_spacetime_type<model_t,model_t::dimension>::signature sp;
		os<<"spacetime type:          "<<abi::__cxa_demangle(typeid(sp).name(),0,0,NULL)<<std::endl;
		os<<std::endl;
		if(model_t::dimension != 0)
		{
		    os<<"continuous helicities:   "<<model_t::continuous_helicities<<std::endl;
		}
		else
		{
		    os<<"continuous helicities:   "<<"NA"<<std::endl;
		}
		os<<std::endl;
		typename get_colour_treatment<model_t>::type ct;
		os<<"colour treatment:        "<<abi::__cxa_demangle(typeid(ct).name(),0,0,NULL)<<std::endl;
		os<<std::endl;
		if(model_t::coloured)
		{
		    os<<"continuous colours:      "<<model_t::continuous_colours<<std::endl;
		}
		else
		{
		    os<<"continuous colours:      "<<"NA"<<std::endl;
		}
		os<<std::endl<<std::endl;
		os<<"particles:"<<std::endl;
		model_wrapper<model_t>::print_particles(os);
		model_wrapper<model_t>::print_families(os);
		os<<"vertices:"<<std::endl;
		model_wrapper<model_t>::print_vertices(os);
		return os;
	    }

	    /// Sets all decay widths to zero.
	    
	    static void discard_widths()
	    {
		width_scheme<model_t>::switch_off();
	    }

	    /// Applies the finite decay widths.

	    static void apply_widths()
	    {
		width_scheme<model_t>::switch_on();
	    }

	    /// Sets the fixed width scheme for all propagators (no width
	    /// included in spacelike propagators, not gauge-invariant).

	    static void set_fixed_width_scheme()
	    {
		width_scheme<model_t>::set_fixed_width_scheme();
	    }

	    /// Sets the complex mass scheme for all propagators (width included
	    /// in spacelike propagators, allows exact gauge-invariance).
	    
	    static void set_complex_mass_scheme()
	    {
		width_scheme<model_t>::set_complex_mass_scheme();
	    }

	    /// Sets the running-width scheme (no complex masses in numerators,
	    /// no widths in spacelike propagators, not gauge-invariant).

	    static void set_running_width_scheme()
	    {
		width_scheme<model_t>::set_running_width_scheme();
	    }

	    /// Switches off spacelike fermion widths.

	    static void set_spacelike_fermion_widths(bool q)
	    {
		fermion_propagator<model_t>::spacelike_width=q;
		anti_fermion_propagator<model_t>::spacelike_width=q;
	    }

	    /// Erases a particle named phi from the static list.
	    /// If no such particle is found, no action is performed.
	    /// Usage may lead to unexpected results as all previously built
	    /// CM_algorithm trees become invalid. For a safer (and faster) way
	    /// to omit particles, use decouple_particle.

	    static void erase_particle(const std::string& phi)
	    {
		model_wrapper<model_t>::erase_particle(phi);
	    }

	    /// Erases the 3-vertices connecting phi1,phi2,phi3.
	    /// If no such vertex is found, no action is performed.
	    /// Usage may lead to unexpected results as all previously built
	    /// CM_algorithm trees become invalid. For a safer (and faster) way
	    /// to omit interactions, use decouple_vertex.

	    static void erase_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3)
	    {
		model_wrapper<model_t>::erase_vertex(phi1,phi2,phi3);
	    }

	    /// Erases the 4-vertices connecting phi1,phi2,str3,str4.
	    /// If no such vertex is found, no action is performed.

	    static void erase_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3,const std::string& phi4)
	    {
		model_wrapper<model_t>::erase_vertex(phi1,phi2,phi3,phi4);
	    }

	    /// Decouples particle named phi from the rest.
	    /// All vertices containing the particle will be disregarded during
	    /// tree evaluations. If phi is not found, no action is performed.

	    static void decouple_particle(const std::string& phi)
	    {
		model_wrapper<model_t>::decouple_particle(phi);
	    }

	    /// Couples particle named phi. 
	    /// If phi is not found or not decoupled, no action is performed.

	    static void couple_particle(const std::string& phi)
	    {
		model_wrapper<model_t>::couple_particle(phi);
	    }

	    /// Decouples the 3-vertices connecting phi1,phi2,phi3.
	    /// If no such vertex is found or no such vertex is coupled, no action is performed.

	    static void decouple_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3)
	    {
		model_wrapper<model_t>::decouple_vertex(phi1,phi2,phi3);
	    }

	    /// Couples the 3-vertices connecting phi1,phi2,phi3.
	    /// If no such vertex is found or no such vertex is decoupled, no action is performed.

	    static void couple_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3)
	    {
		model_wrapper<model_t>::couple_vertex(phi1,phi2,phi3);
	    }

	    /// Decouples the 4-vertices connecting phi1,phi2,phi3,phi4.

	    static void decouple_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3,const std::string& phi4)
	    {
		model_wrapper<model_t>::decouple_vertex(phi1,phi2,phi3,phi4);
	    }

	    /// Couples the 4-vertices connecting phi1,phi2,phi3,phi4.

	    static void couple_vertex(const std::string& phi1,const std::string& phi2,const std::string& phi3,const std::string& phi4)
	    {
		model_wrapper<model_t>::couple_vertex(phi1,phi2,phi3,phi4);
	    }

	    /// decouples all vertices with zero couplings and couples all vertices with nonzero couplings

	    static void decouple_zero_vertices()
	    {
		model_wrapper<model_t>::decouple_zero_vertices();
	    }

	    /// Replaces the propagator of particle named 'phi' by prop_t.

	    template<template<class M>class prop_t>static void set_propagator(const std::string& phi)
	    {
		model_wrapper<model_t>::template set_propagator< prop_t<model_t> >(phi);
	    }

	    /// Outputs the model's fusion rules used in the algorithm.

	    static std::ostream& print_fusion_rules(std::ostream& os)
	    {
		return model_wrapper<model_t>::print_fusion_rules(os);
	    }

	    /// Assigns mass *m to the particle named phi.

	    static void set_mass(const std::string& phi,void* m)
	    {
		model_wrapper<model_t>::set_mass(phi,static_cast<const typename model_t::value_type*>(m));
	    }

	    /// Makes phi massless.

	    static void set_massless(const std::string& phi)
	    {
		model_wrapper<model_t>::set_massless(phi);
	    }

	    /// Assigns width *w to the particle named phi.

	    static void set_width(const std::string& phi,void* w)
	    {
		model_wrapper<model_t>::set_width(phi,static_cast<const typename model_t::value_type*>(w));
	    }

	    /// Makes phi widthless.

	    static void set_widthless(const std::string& phi)
	    {
		model_wrapper<model_t>::set_widthless(phi);
	    }

	    /// Sets the QCD scale (empty when not overloaded).

	    static void set_QCD_scale(){}

	    /// Sets the strong coupling (empty when not overloaded).

	    static void set_alpha_s(){}

	protected:

	    /* Protected constructor, is called when a specific model is created
	     * and opens and initialises the model log file: */

	    model()
	    {
		/* Count the created model object: */

		++obj_counter;
		
		/* Retrieve the name of the subclassed model: */
		
		modelname=abi::__cxa_demangle(typeid(instance).name(),0,0,NULL);
		modelname.erase(modelname.end()-1);

		/* Attempt to open a logfile: */

		std::string barename=modelname;
		std::string::reverse_iterator it=barename.rbegin();
		while(it!=barename.rend() and *it!=':')
		{
		    ++it;
		}
		barename.erase(0,barename.rend()-it);
		fs.open(barename.append(".log").c_str());
		if(fs.is_open())
		{
		    /* General information writing: */
		    
		    fs<<"Log file for model "<<modelname<<", created on "<<__DATE__<<", at "<<__TIME__<<std::endl;
		    fs<<std::endl;

		    /* Writing the compile-time model information: */

		    fs<<"dimension:               "<<model_t::dimension<<std::endl;
		    fs<<std::endl;
		    typename model_t::value_type x;
		    fs<<"value type:              "<<abi::__cxa_demangle(typeid(x).name(),0,0,NULL)<<std::endl;
		    fs<<std::endl;
		    typename get_spacetime_type<model_t,model_t::dimension>::signature sp;
		    fs<<"spacetime type:          "<<abi::__cxa_demangle(typeid(sp).name(),0,0,NULL)<<std::endl;
		    fs<<std::endl;
		    if(model_t::dimension != 0)
		    {
			fs<<"continuous helicities:   "<<model_t::continuous_helicities<<std::endl;
		    }
		    else
		    {
			fs<<"continuous helicities:   "<<"NA"<<std::endl;
		    }
		    fs<<std::endl;
		    typename get_colour_treatment<model_t,model_t::coloured>::type ct;
		    fs<<"colour treatment:        "<<abi::__cxa_demangle(typeid(ct).name(),0,0,NULL)<<std::endl;
		    fs<<std::endl;
		    if(!model_t::coloured)
		    {
			fs<<"continuous colours:      "<<get_colour_treatment<model_t>::continuous_colours<<std::endl;
		    }
		    else
		    {
			fs<<"continuous colours:      "<<"NA"<<std::endl;
		    }
		    fs<<std::endl<<std::endl;
		}
		else
		{
		    /* Output error message: */

		    std::cout<<"Failed to open log file for "<<modelname<<std::endl;
		    std::cout<<"No logging will be performed."<<std::endl;
		}
	    }

	private:

	    /* Singleton instance: */

	    static model_t* instance;
	    
	    /* Log stream: */
	    
	    static std::ofstream fs;
	    
	    /* Subclassed model name: */
	    
	    static std::string modelname;

	    /* Object counter: */

	    static int obj_counter;
    };
    template<class model_t> std::ofstream model<model_t>::fs;
    template<class model_t> model_t* model<model_t>::instance=NULL;
    template<class model_t>std::string model<model_t>::modelname;
    template<class model_t>int model<model_t>::obj_counter=0;
}

#endif /*CAMGEN_MODEL_H_*/

