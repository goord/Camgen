//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_GGG_H_
#define CAMGEN_GGG_H_

namespace Camgen
{
    /* Specialisation of the Yang-Mills 3-vertex composed with the
     * SU(N)-structure constants (e.g. the 3-gluon vertex in N-colour QCD): */

    template<std::size_t N,class model_t>class evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);
	    
	    /* SU(N)-group type definition: */
	    
	    DEFINE_GROUP_TYPE(model_t,SU<N>);

	    /* Convenient integer constants: */

	    static const std::size_t vector_size=model_t::dimension;
	    static const std::size_t gluon_size=vector_size*(N*N-1);

	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
		group_type::initialise();
	    }

	    /* Utility function returning the index ranges of the interacting
	     * subamplitudes: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n<3)
		{
		    v.push_back(N*N-1);
		    v.push_back(vector_size);
		}
		return v;
	    }

	    /* Recursive relations: */

	    static void first(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_0.range()<(N*N-1)*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<(N*N-1)*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<(N*N-1)*vector_size),"iterator 2 out of range");

		/* Storing inner products in memory: */

		value_type nul(0,0);
		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<a;++b)
		    {
			c3[a][b]+=spacetime_type::dot(A_1,A_2);
			A_2+=vector_size;
		    }
		    c1[a]=(r_value_type)2*spacetime_type::dot(P_1,A_2)+spacetime_type::dot(P_2,A_2);
		    c2[a]=(r_value_type)2*spacetime_type::dot(P_2,A_1)+spacetime_type::dot(P_1,A_1);
		    A_2+=vector_size;
		    for(size_type b=a+1;b<N*N-1;++b)
		    {
			c3[b][a]=-spacetime_type::dot(A_1,A_2);
			A_2+=vector_size;
		    }
		    A_2-=gluon_size;
		    A_1+=vector_size;
		}
		A_1-=gluon_size;
		iterator B_1=A_1;
		iterator B_2=A_2;

		/* Performing contractions: */

		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<N*N-1;++b)
		    {
			for(size_type c=0;c<b;++c)
			{
			    if(group_type::structure_constant(a,b,c)!=nul)
			    {
				for(size_type mu=0;mu<vector_size;++mu)
				{
				    A_0[mu]+=factor*C_0*(group_type::structure_constant(a,b,c))*(-c1[c]*A_1[mu]+c1[b]*B_1[mu]+c3[b][c]*(P_1[mu]-P_2[mu])+c2[b]*A_2[mu]-c2[c]*B_2[mu]);
				}
			    }
			    A_2+=vector_size;
			    B_1+=vector_size;
			}
			A_1+=vector_size;
			B_2+=vector_size;
			A_2-=b*vector_size;
			B_1-=b*vector_size;
		    }
		    A_0+=vector_size;
		    A_1-=gluon_size;
		    B_2-=gluon_size;
		}
		A_0-=gluon_size;

	    }
	    static void second(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_0.range()<(N*N-1)*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<(N*N-1)*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<(N*N-1)*vector_size),"iterator 2 out of range");

		/* Storing inner products in memory: */

		value_type nul(0,0);
		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<a;++b)
		    {
			c3[a][b]+=spacetime_type::dot(A_0,A_2);
			A_2+=vector_size;
		    }
		    c1[a]=(r_value_type)2*spacetime_type::dot(P_0,A_2)+spacetime_type::dot(P_2,A_2);
		    c2[a]=(r_value_type)2*spacetime_type::dot(P_2,A_0)+spacetime_type::dot(P_0,A_0);
		    A_2+=vector_size;
		    for(size_type b=a+1;b<N*N-1;++b)
		    {
			c3[b][a]=-spacetime_type::dot(A_0,A_2);
			A_2+=vector_size;
		    }
		    A_2-=gluon_size;
		    A_0+=vector_size;
		}
		A_0-=gluon_size;
		iterator B_0=A_0;
		iterator B_2=A_2;

		/* Performing contractions: */

		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<N*N-1;++b)
		    {
			for(size_type c=0;c<b;++c)
			{
			    if(group_type::structure_constant(a,b,c)!=nul)
			    {
				for(size_type mu=0;mu<vector_size;++mu)
				{
				    A_1[mu]+=factor*C_0*(group_type::structure_constant(a,b,c))*(-c1[c]*A_0[mu]+c1[b]*B_0[mu]+c3[b][c]*(P_0[mu]-P_2[mu])+c2[b]*A_2[mu]-c2[c]*B_2[mu]);
				}
			    }
			    A_2+=vector_size;
			    B_0+=vector_size;
			}
			A_0+=vector_size;
			B_2+=vector_size;
			A_2-=b*vector_size;
			B_0-=b*vector_size;
		    }
		    A_1+=vector_size;
		    A_0-=gluon_size;
		    B_2-=gluon_size;
		}
		A_1-=gluon_size;

	    }
	    static void third(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_0.range()<(N*N-1)*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<(N*N-1)*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<(N*N-1)*vector_size),"iterator 2 out of range");

		/* Storing inner products in memory: */

		value_type nul(0,0);
		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<a;++b)
		    {
			c3[a][b]+=spacetime_type::dot(A_0,A_1);
			A_1+=vector_size;
		    }
		    c1[a]=(r_value_type)2*spacetime_type::dot(P_0,A_1)+spacetime_type::dot(P_1,A_1);
		    c2[a]=(r_value_type)2*spacetime_type::dot(P_1,A_0)+spacetime_type::dot(P_0,A_0);
		    A_1+=vector_size;
		    for(size_type b=a+1;b<N*N-1;++b)
		    {
			c3[b][a]=-spacetime_type::dot(A_0,A_1);
			A_1+=vector_size;
		    }
		    A_1-=gluon_size;
		    A_0+=vector_size;
		}
		A_0-=gluon_size;
		iterator B_0=A_0;
		iterator B_1=A_1;

		/* Performing contractions: */

		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<N*N-1;++b)
		    {
			for(size_type c=0;c<b;++c)
			{
			    if(group_type::structure_constant(a,b,c)!=nul)
			    {
				for(size_type mu=0;mu<vector_size;++mu)
				{
				    A_2[mu]+=factor*C_0*(group_type::structure_constant(a,b,c))*(-c1[c]*A_0[mu]+c1[b]*B_0[mu]+c3[b][c]*(P_0[mu]-P_1[mu])+c2[b]*A_1[mu]-c2[c]*B_1[mu]);
				}
			    }
			    A_1+=vector_size;
			    B_0+=vector_size;
			}
			A_0+=vector_size;
			B_1+=vector_size;
			A_1-=b*vector_size;
			B_0-=b*vector_size;
		    }
		    A_2+=vector_size;
		    A_0-=gluon_size;
		    B_1-=gluon_size;
		}
		A_2-=gluon_size;

	    }
	    static void fourth(ARG_LIST){}
	private:
	    static value_type c1[N*N-1];
	    static value_type c2[N*N-1];
	    static value_type c3[N*N-1][N*N-1];
    };
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::vector_size;
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::gluon_size;
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::c1[N*N-1]={0};
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::c2[N*N-1]={0};
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::f<SU<N>,0,1,2>,vvv<model_t> > >::c3[N*N-1][N*N-1]={{0}};

    /* Specialisation evaluate class template for the gluon 3-vertex with the
     * gluons in the colour-flow representation: */

    template<std::size_t N,class model_t>class evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);
	    
	    /* Group type definition: */
	    
	    DEFINE_GROUP_TYPE(model_t,SU<N>);

	    /* Useful constants: */

	    static const std::size_t vector_size=model_t::dimension;
	    static const std::size_t quark_size=N*vector_size;
	    static const std::size_t gluon_size=N*quark_size;

	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
	    }

	    /* Checking function returning the index ranges of interacting
	     * subamplitude tensors: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n<3)
		{
		    v.push_back(N);
		    v.push_back(N);
		    v.push_back(vector_size);
		}
		return v;
	    }

	    /* Recursive relations: */

	    static void first(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_0.range()<N*N*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<N*N*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<N*N*vector_size),"iterator 2 out of range");

		/* Storing inner product combinations in memory: */

		for(size_type a=0;a<N;++a)
		{
		    for(size_type b=0;b<N;++b)
		    {
			c1[a][b]=(r_value_type)2*spacetime_type::dot(P_1,A_2)+spacetime_type::dot(P_2,A_2);
			c2[a][b]=(r_value_type)2*spacetime_type::dot(P_2,A_1)+spacetime_type::dot(P_1,A_1);	
			A_1+=vector_size;
			A_2+=vector_size;
		    }
		}
		A_1-=gluon_size;
		A_2-=gluon_size;
		
		/* Performing contractions: */
		
		iterator B_1=A_1;
		iterator B_2=A_2;
		value_type z(0,-1);
		z*=(factor*C_0);
		value_type temp;
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    temp=spacetime_type::dot(A_1,B_2)-spacetime_type::dot(A_2,B_1);
			    for(size_type mu=0;mu<vector_size;++mu)
			    {
				A_0[mu]+=z*(c1[a][c]*B_1[mu]-A_1[mu]*c1[c][b]+c2[a][c]*B_2[mu]-A_2[mu]*c2[c][b]+temp*(P_1[mu]-P_2[mu]));
			    }
			    A_0+=vector_size;
			    B_1+=vector_size;
			    B_2+=vector_size;
			}
			A_0-=quark_size;
			A_1+=vector_size;
			A_2+=vector_size;
		    }
		    A_0+=quark_size;
		    B_1-=gluon_size;
		    B_2-=gluon_size;
		}
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_2-=gluon_size;

	    }
	    static void second(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_0.range()<N*N*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<N*N*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<N*N*vector_size),"iterator 2 out of range");

		/* Storing inner product combinations in memory: */

		for(size_type a=0;a<N;++a)
		{
		    for(size_type b=0;b<N;++b)
		    {
			c1[a][b]=(r_value_type)2*spacetime_type::dot(P_0,A_2)+spacetime_type::dot(P_2,A_2);
			c2[a][b]=(r_value_type)2*spacetime_type::dot(P_2,A_0)+spacetime_type::dot(P_0,A_0);	
			A_0+=vector_size;
			A_2+=vector_size;
		    }
		}
		A_0-=gluon_size;
		A_2-=gluon_size;
		iterator B_0=A_0;
		iterator B_2=A_2;
		value_type z(0,-1);
		z*=(factor*C_0);
		value_type temp;
		
		/* Performing contractions: */
		
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    temp=spacetime_type::dot(A_0,B_2)-spacetime_type::dot(A_2,B_0);
			    for(size_type mu=0;mu<vector_size;++mu)
			    {
				A_1[mu]+=z*(c1[a][c]*B_0[mu]-A_0[mu]*c1[c][b]+c2[a][c]*B_2[mu]-A_2[mu]*c2[c][b]+temp*(P_0[mu]-P_2[mu]));
			    }
			    A_1+=vector_size;
			    B_0+=vector_size;
			    B_2+=vector_size;
			}
			A_1-=quark_size;
			A_0+=vector_size;
			A_2+=vector_size;
		    }
		    A_1+=quark_size;
		    B_0-=gluon_size;
		    B_2-=gluon_size;
		}
		A_1-=gluon_size;
		A_0-=gluon_size;
		A_2-=gluon_size;

	    }
	    static void third(ARG_LIST)
	    {

		CAMGEN_ERROR_IF((A_0.range()<N*N*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<N*N*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<N*N*vector_size),"iterator 2 out of range");

		/* Storing inner product combinations in memory: */

		for(size_type a=0;a<N;++a)
		{
		    for(size_type b=0;b<N;++b)
		    {
			c1[a][b]=(r_value_type)2*spacetime_type::dot(P_0,A_1)+spacetime_type::dot(P_1,A_1);
			c2[a][b]=(r_value_type)2*spacetime_type::dot(P_1,A_0)+spacetime_type::dot(P_0,A_0);	
			A_0+=vector_size;
			A_1+=vector_size;
		    }
		}
		A_0-=gluon_size;
		A_1-=gluon_size;
		iterator B_0=A_0;
		iterator B_1=A_1;
		value_type z(0,-1);
		z*=(factor*C_0);
		value_type temp;
		
		/* Performing contractions: */
		
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    temp=spacetime_type::dot(A_0,B_1)-spacetime_type::dot(A_1,B_0);
			    for(size_type mu=0;mu<vector_size;++mu)
			    {
				A_2[mu]+=z*(c1[a][c]*B_0[mu]-A_0[mu]*c1[c][b]+c2[a][c]*B_1[mu]-A_1[mu]*c2[c][b]+temp*(P_0[mu]-P_1[mu]));
			    }
			    A_2+=vector_size;
			    B_0+=vector_size;
			    B_1+=vector_size;
			}
			A_2-=quark_size;
			A_0+=vector_size;
			A_1+=vector_size;
		    }
		    A_2+=quark_size;
		    B_0-=gluon_size;
		    B_1-=gluon_size;
		}
		A_2-=gluon_size;
		A_0-=gluon_size;
		A_1-=gluon_size;

	    }
	    static void fourth(ARG_LIST){}
	private:
	    static value_type c1[N][N];
	    static value_type c2[N][N];
    };
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >::vector_size;
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >::quark_size;
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >::gluon_size;
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >::c1[N][N]={{0}};
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::CF_f<N,0,1,2>,vvv<model_t> > >::c2[N][N]={{0}};
}

#endif /*CAMGEN_GGG_H_*/
