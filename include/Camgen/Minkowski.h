//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_MINKOWSKI_H_
#define CAMGEN_MINKOWSKI_H_

#include <Camgen/spacetime.h>
#include <Camgen/vector.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
 * Minkowski spacetime definition. For models in dimensions > 0, a spacetime type    *
 * definition is required for the inner products of momenta and Lorents vectors. A   *
 * spacetime class wraps an implementation metafunction, taking the value type and   *
 * dimension as parameters and defining the inner products between real- or          *
 * complex-values Lorentz vectors. Moreover, the spacetime implementation class      *
 * contains functions for constructing the base vectors k_0 and k_1 for massles basic*
 * spinors and functions that construct spin vectors for massive spinors. The        *
 * Minkowski metric in Camgen has signature (+ - - -). 			     	     *
 *                               						     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    class Minkowski_type
    {
	public:
	    
	    /* Implementation metafunction class of Minkowksi spacetime. This is
	     * a class derived from the spacetime base class template: */

	    template<class value_t,std::size_t dim>class implementation: public spacetime<value_t,implementation<value_t,dim>,dim>
	    {
	    	public:
		    
		    /* Type definition referring to the base class: */
		    
		    typedef spacetime<value_t,implementation<value_t,dim>,dim> base_type;
		    
		    /* Static constant boolean denoting whether the metric is
		     * diagonal (for optimisation purposes): */

		    static const bool diagonal=true;

		    /* Time (energy) component definition: */

		    static const std::size_t timelike_direction=0;
		    
		    /* Implementation of the metric tensor (entries are by
		     * default zero): */

		    static void fill_metric()
		    {
			base_type::g[0][0]=(value_t)(1);
			for(std::size_t i=1;i<dim;++i)
			{
			    base_type::g[i][i]=(value_t)(-1);
			}
		    }

		    /* Separate implementation of the inner product between
		     * vectors. The type member of product_type<T1,T2> assigns
		     * the correct value type (real if T1 and T2 are both
		     * real,complex if one of them is): */

		    template<class vector_t1,class vector_t2>static typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type dot(const vector_t1& t1,const vector_t2& t2)
		    {
			typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type result=t1[0]*t2[0];
			for(std::size_t i=1;i<dim;++i)
			{
			    result-=t1[i]*t2[i];
			}
			return result;
		    }

		    /* Separate implementation of the space-like product between
		     * vectors. The type member of product_type<T1,T2> assigns
		     * the correct value type (real if T1 and T2 are both
		     * real,complex if one of them is): */

		    template<class vector_t1,class vector_t2>static typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type space_dot(const vector_t1& t1,const vector_t2& t2)
		    {
			typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type result=t1[1]*t2[1];
			for(std::size_t i=2;i<dim;++i)
			{
			    result+=t1[i]*t2[i];
			}
			return result;
		    }
		  
		    /* Function defining base vectors for the construction of
		     * the basic massless spinors, w.r.t. the beam direction: */

		    template<int beam_dir>static void make_base_vectors(vector<value_t,dim>& k_0,vector<value_t,dim>& k_1,vector<value_t,dim>& k_2)
		    {
			CAMGEN_ERROR_IF((beam_dir<-dim+1 or beam_dir>dim-1 or beam_dir==0),"invalid beam direction chosen");

			k_0.assign((value_t)0);
			k_1.assign((value_t)0);
			k_2.assign((value_t)0);
			
			if(beam_dir>2 or beam_dir<-2)
			{
			    k_0[0]=(value_t)1;
			    k_0[1]=(value_t)1;
			    k_1[2]=(value_t)1;
			    k_2[0]=(value_t)1;
			    k_2[2]=(value_t)1;
			    return;
			}
			if(beam_dir==1 or beam_dir==-1)
			{
			    k_0[0]=(value_t)1;
			    k_0[2]=(value_t)1;
			    k_1[3]=(value_t)1;
			    k_2[0]=(value_t)1;
			    k_2[3]=(value_t)1;
			    return;
			}
			else if(beam_dir==2 or beam_dir==-2)
			{
			    k_0[0]=(value_t)1;
			    k_0[1]=(value_t)1;
			    k_1[3]=(value_t)1;
			    k_2[0]=(value_t)1;
			    k_2[3]=(value_t)1;
			    return;	
			}
		    }
		    
		    /* Function defining the construction of the helicity spin vector for
		     * massive spinors. This spin vector is the longitudinal
		     * polarisation vector to the momentum, and will hence
		     * produce NAN in the rest frame of the fermion: */

		    static vector<value_t,dim>& make_helicity_spinvec(vector<value_t,dim>& s,const vector<value_t,dim>& p,const value_t& m)
		    {
			CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
			CAMGEN_MESSAGE_IF((p[0]==m or p[0]==-m),"center-of-mass momentum encountered");

			value_t P=std::sqrt(p[0]*p[0]-m*m);
			value_t t=p[0]/(P*m);
			s[0]=P/m;
			for(std::size_t i=1;i<dim;++i)
			{
			    s[i]=t*p[i];
			}
			return s;
		    }

		    /* Function defining the construction of the Kleiss-Stirling
		     * spin vector. It is defined as
		     *
		     *			s=(1/m)*p - (m/(k_0.p))*k_0
		     * 
		     * where k_0 is the first vector defined by the
		     * make_base_vectors() routine. */

		    template<int beam_dir>static vector<value_t,dim>& make_KS_spinvec(vector<value_t,dim>& s,const vector<value_t,dim>& p,const value_t& m)
		    {
			int n=1;
			if(beam_dir==1 or beam_dir==-1)
			{
			    n=2;
			}

			CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
			CAMGEN_MESSAGE_IF((p[0]==p[n]),"zero-mass momentum parallel to k_0 encountered");
			
			s=p/m;
			value_t f=m/(p[0]-p[n]);
			s[0]-=f;
			s[n]-=f;
			return s;
		    }

		    /* Function defining the construction of the polarised spin
		     * vector. It is defined as
		     *
		     *              s=e - (e.p/(m*(E + m))*(p + p0)
		     *
		     * where e is the spacelike unit vector in the polarisation
		     * direction and p0 the rest-frame four-momentum. Note that
		     * in the rest-frame s=e.
		     * */

		    template<int beam_dir>static vector<value_t,dim>& make_polarised_spinvec(vector<value_t,dim>& s,const vector<value_t,dim>& p,const value_t& m)
		    {
			CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
			CAMGEN_MESSAGE_IF((p[0]==m),"center-of-mass momentum encountered");
			
			if(beam_dir>0)
			{
			    s[0]=p[beam_dir]/m;
			    value_t a=s[0]/(p[0]+m);
			    for(std::size_t i=1;i<dim;++i)
			    {
				s[i]=a*p[i];
			    }
			    s[beam_dir]+=((value_t)1);
			}
			if(beam_dir<0)
			{
			    s[0]=-p[-beam_dir]/m;
			    value_t a=s[0]/(p[0]+m);
			    for(std::size_t i=1;i<dim;++i)
			    {
				s[i]=a*p[i];
			    }
			    s[-beam_dir]-=((value_t)1);
			}
			return s;
		    }

		    /* Momentum-splitting function for the construction of
		     * massive polarisation vector: */

		    static void split_momentum(vector<value_t,dim>& p1,vector<value_t,dim>& p2,const vector<value_t,dim>& p,const value_t& m)
		    {
			CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
			CAMGEN_MESSAGE_IF((p[0]==m or p[0]==-m),"center-of-mass momentum encountered");
			
			value_t msq=m*m;
			value_t pvec=std::sqrt(p[0]*p[0]-msq);
			value_t a=0.5*msq/(pvec*(p[0]-pvec));
			p1[0]=a*pvec;
			for(std::size_t i=1;i<dim;++i)
			{
			    p1[i]=a*p[i];
			}
			p2=p;
			p2-=p1;
		    }
	};
    };
    template<class value_t,std::size_t dim>const std::size_t Minkowski_type::implementation<value_t,dim>::timelike_direction;
    template<class value_t,std::size_t dim>const bool Minkowski_type::implementation<value_t,dim>::diagonal;

    /* For optimisation purposes, a separate specialisation of the
     * four-dimensional case, unrolling the loops: */

    template<class value_t>class Minkowski_type::implementation<value_t,4>: public spacetime<value_t,Minkowski_type::implementation<value_t,4>,4>
    {
	public:		    
	    /* Type definition referring to the base class: */
		   
	    typedef spacetime<value_t,implementation<value_t,4>,4> base_type;		    
	    
	    /* Static constant boolean denoting whether the metric is
	     * diagonal (for optimisation purposes): */

	    static const bool diagonal=true;

	    /* Time (energy) component definition: */

	    static const std::size_t timelike_direction=0;
	
	    /* Implementation of the metric tensor (entries are by
	     * default zero): */

	    static void fill_metric()
	    {
		base_type::g[0][0]=(value_t)1;
		base_type::g[1][1]=-(value_t)1;
		base_type::g[2][2]=-(value_t)1;
		base_type::g[3][3]=-(value_t)1;
	    }	
	    
	    /* Separate implementation of the inner product between vectors. The
	     * type member of product_type<T1,T2> assigns the correct value type
	     * (real if T1 and T2 are both real,complex if one of them is): */

	    template<class vector_t1,class vector_t2>static typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type dot(const vector_t1& t1,const vector_t2& t2)
	    {
		return (t1[0]*t2[0]-t1[1]*t2[1]-t1[2]*t2[2]-t1[3]*t2[3]);
	    }
	    
	    /* Separate implementation of the space-like inner product between
	     * vectors. The type member of product_type<T1,T2> assigns the
	     * correct value type (real if T1 and T2 are both real,complex if
	     * one of them is): */

	    template<class vector_t1,class vector_t2>static typename product_type<value_t,typename vector_t1::value_type,typename vector_t2::value_type>::value_type space_dot(const vector_t1& t1,const vector_t2& t2)
	    {
		return (t1[1]*t2[1]+t1[2]*t2[2]+t1[3]*t2[3]);
	    }
	    
	    /* Function defining base vectors for the construction of
	     * the basic massless spinors, w.r.t. the beam direction: */

	    template<int beam_dir>static void make_base_vectors(vector<value_t,4>& k_0,vector<value_t,4>& k_1,vector<value_t,4>& k_2)
	    {
		CAMGEN_ERROR_IF((beam_dir<-3 or beam_dir>3 or beam_dir==0),"invalid beam direction chosen");

		k_0.assign((value_t)0);
		k_1.assign((value_t)0);
		k_2.assign((value_t)0);
		if(beam_dir==3 or beam_dir==-3)
		{
		    k_0[0]=(value_t)1;
		    k_0[1]=(value_t)1;
		    k_1[2]=(value_t)1;
		    k_2[0]=(value_t)1;
		    k_2[2]=(value_t)1;
		    return;
		}
		if(beam_dir==1 or beam_dir==-1)
		{
		    k_0[0]=(value_t)1;
		    k_0[2]=(value_t)1;
		    k_1[3]=(value_t)1;
		    k_2[0]=(value_t)1;
		    k_2[3]=(value_t)1;
		    return;
		}
		else if(beam_dir==2 or beam_dir==-2)
		{
		    k_0[0]=(value_t)1;
		    k_0[1]=(value_t)1;
		    k_1[3]=(value_t)1;
		    k_2[0]=(value_t)1;
		    k_2[3]=(value_t)1;
		    return;	
		}
	    }		    
	    
	    /* Function defining the construction of the helicity spin vector for
	     * massive spinors. This spin vector is the longitudinal
             * polarisation vector to the momentum, and will hence
             * produce NAN in the rest frame of the fermion: */

	    static vector<value_t,4>& make_helicity_spinvec(vector<value_t,4>& s,const vector<value_t,4>& p,const value_t& m)
	    {
		CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
		CAMGEN_MESSAGE_IF((p[0]==m or p[0]==-m),"center-of-mass momentum encountered");

		value_t P=std::sqrt(p[0]*p[0]-m*m);
		value_t t=p[0]/(P*m);
		s[0]=P/m;
		s[1]=t*p[1];
		s[2]=t*p[2];
		s[3]=t*p[3];
		return s;
	    }

	    /* Function defining the construction of the Kleiss-Stirling
	     * spin vector. It is defined as
	     *
	     *			s=(1/m)*p - (m/(k_0.p))*k_0
	     * 
	     * where k_0 is the first vector defined by the
	     * make_base_vectors() routine. */

	    template<int beam_dir>static vector<value_t,4>& make_KS_spinvec(vector<value_t,4>& s,const vector<value_t,4>& p,const value_t& m)
	    {
		int n=1;
		if(beam_dir==1 or beam_dir==-1)
		{
		    n=2;
		}

		CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
		CAMGEN_MESSAGE_IF((p[0]==p[n]),"zero-mass momentum parallel to k_0 encountered");

		s=p;
		s/=m;
		value_t f=m/(p[0]-p[n]);
		s[0]-=f;
		s[n]-=f;
		return s;
	    }

	    /* Function defining the construction of the polarised spin
	     * vector. It is defined as
	     *
	     *              s=e - (e.p/(m*(E + m))*(p + p0)
	     *
	     * where e is the spacelike unit vector in the polarisation
	     * direction and p0 the rest-frame four-momentum. Note that
	     * in the rest-frame s=e.
	     * */

	    template<int beam_dir>static vector<value_t,4>& make_polarised_spinvec(vector<value_t,4>& s,const vector<value_t,4>& p,const value_t& m)
	    {
		CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
		CAMGEN_MESSAGE_IF((p[0]==-m),"center-of-mass momentum encountered");
			
		if(beam_dir>0)
		{
		    s[0]=p[beam_dir]/m;
		    value_t a=s[0]/(p[0]+m);
		    s[1]=a*p[1];
		    s[2]=a*p[2];
		    s[3]=a*p[3];
		    s[beam_dir]+=((value_t)1);
		}
		if(beam_dir<0)
		{
		    s[0]=-p[-beam_dir]/m;
		    value_t a=s[0]/(p[0]+m);
		    s[1]=a*p[1];
		    s[2]=a*p[2];
		    s[3]=a*p[3];
		    s[-beam_dir]-=((value_t)1);
		}
		return s;
	    }

	    /* Momentum-splitting function for the construction of
	     * massive polarisation vector: */

	    static void split_momentum(vector<value_t,4>& p1,vector<value_t,4>& p2,const vector<value_t,4>& p,const value_t& m)
	    {
		CAMGEN_MESSAGE_IF((m<=(value_t)0),"zero or negative mass encountered");
		CAMGEN_MESSAGE_IF((p[0]==m or p[0]==-m),"center-of-mass momentum encountered");
			
		value_t msq=m*m;
		value_t pvec=std::sqrt(p[0]*p[0]-msq);
		value_t a=(value_t)0.5*msq/(pvec*(p[0]-pvec));
		p1[0]=a*pvec;
		p1[1]=a*p[1];
		p1[2]=a*p[2];
		p1[3]=a*p[3];
		p2=p;
		p2-=p1;
	    }
    };
    template<class value_t>const std::size_t Minkowski_type::implementation<value_t,4>::timelike_direction;
    template<class value_t>const bool Minkowski_type::implementation<value_t,4>::diagonal;
}

#endif /*CAMGEN_MINKOWSKI_H_*/


