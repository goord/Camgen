//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_GGGG_H_
#define CAMGEN_GGGG_H_

namespace Camgen
{
    /* Declaration and definition of the four-gluon vertex with the gluons in
     * the adjoint representation: */

    template<std::size_t N,class model_t>class compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >
    {
	private:

	    /* Gluon size defnition: */

	    static const std::size_t gluon_size=(model_t::dimension)*(N*N-1);

	public:

	    /* Rank of the vertex: */

	    static const std::size_t rank=4;

	    /* Number of couplings used by the vertex: */

	    static const std::size_t params=1;
	    
	    /* Boolean denoting whether the Feynman rule depends on momenta: */
	    
	    static const bool p_dependent=false;

	    /* Boolean denoting whether the Feynman rule is fermionic: */

	    static const bool fermionic=false;
	    
	    /* Vertex tensor size definition: */
	    
	    static const std::size_t tensor_size=gluon_size*gluon_size*gluon_size*gluon_size;

	    /* Subamplitude tensor sizes: */

	    static const std::size_t sizes[4];

	    /* Function printing the formula of the Feynman rule: */

	    static const std::string formula;
    };
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >::rank;
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >::params;
    template<std::size_t N,class model_t>const bool compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >::p_dependent;
    template<std::size_t N,class model_t>const bool compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >::fermionic;
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >::tensor_size;
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >::sizes[4]={(N*N-1)*model_t::dimension,(N*N-1)*model_t::dimension,(N*N-1)*model_t::dimension,(N*N-1)*model_t::dimension};
    template<std::size_t N,class model_t>const std::string compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> >::formula="f(a1,a2,e)f(a3,a4,e)(g(mu1,mu3)g(mu2,mu4)-g(mu1,mu4)g(mu2,mu3))+2<->3+2<->4";

    namespace colour_tensor
    {
	/* Specialisation of the ff_contr helper class in the case of the
	 * 4-gluon vertex in the adjoint representation: */

	template<std::size_t N,class model_t>class ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>
	{
	    public:

		/* Some usefule type definitions: */

		typedef typename evaluate< vvvv<model_t> >::value_type value_type;
		typedef typename evaluate< vvvv<model_t> >::size_type size_type;
		
		/* Group type definition, source of the structure constants: */
		
		typedef typename SU<N>::template implementation<typename model_t::value_type> group_type;

	    private:

		/* Useful constants: */

		static const std::size_t group_rank=N*N-1;
		static const std::size_t dim=model_t::dimension;
	    
	    public:

		/* Initialisation phase: */

		static void initialise()
		{
		    if(!initialised)
		    {
			/* Filling the indices array with vectors of contracted
			 * structure constants: */

			group_type::initialise();
			value_type nul=0.;
			value_type f;
			for(size_type a=0;a<group_rank;++a)
			{
			    for(size_type b=0;b<=a;++b)
			    {
				for(size_type c=0;c<group_rank;++c)
				{
				    for(size_type d=0;d<=c;++d)
				    {
					f=nul;
					for(size_type e=0;e<group_rank;++e)
					{
					    f+=((group_type::structure_constant(e,a,c))*(group_type::structure_constant(e,b,d)));
					    f+=((group_type::structure_constant(e,a,d))*(group_type::structure_constant(e,b,c)));
					}
					if(f!=nul)
					{
					    values.push_back(f);
					    indices[0].push_back(a);
					    indices[1].push_back(b);
					    indices[2].push_back(c);
					    indices[3].push_back(d);
					}
				    }
				}
			    }
			}

			/* Creating the jumps between indices: */

			for(size_type i=0;i<4;++i)
			{
			    jumps[i].resize(size());
			    jumps[i][0]=dim*indices[i][0];
			    for(size_type j=1;j<size();++j)
			    {
				jumps[i][j]=dim*(indices[i][j]-indices[i][j-1]);
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

		/* Nonzero colour structure values: */

		static value_type value(size_type i)
		{
		    return values[i];
		}

		/* Indices of the nonzero colour structure components: */

		static size_type index(size_type i,size_type j)
		{
		    return indices[i][j];
		}

		/* Offsets to nonzero colour structure components: */

		static int jump(size_type i,size_type j)
		{
		    return jumps[i][j];
		}

		/* Offsets to the first nonzero colour structure components: */

		static size_type reset(size_type i)
		{
		    return dim*(indices[i].back());
		}

		/* Output function: */

		static void print()
		{
		    initialise();
		    for(size_type i=0;i<values.size();++i)
		    {
			std::cout<<indices[0][i]<<","<<indices[1][i]<<","<<indices[2][i]<<","<<indices[3][i]<<"\t\t"<<values[i]<<std::endl;
		    }
		}

	    private:

		/* Values of nonzero components: */

		static std::vector<value_type> values;
		
		/* Indices of nonzero components: */
		
		static std::vector<size_type>indices[4];

		/* Linear offsets to nonzero components: */

		static std::vector<int>jumps[4];
		
		/* Initialisation tag: */
		
		static bool initialised;
	};
	template<std::size_t N,class model_t>const std::size_t ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::group_rank;
	template<std::size_t N,class model_t>const std::size_t ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::dim;
	template<std::size_t N,class model_t>bool ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::initialised=false;
	template<std::size_t N,class model_t>std::vector<typename ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::value_type> ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::values;
	template<std::size_t N,class model_t>std::vector<typename ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::size_type> ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::indices[4];
	template<std::size_t N,class model_t>std::vector<int>ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3>::jumps[4];
    }

    /* Evaluate class template specialisation for the four-gluon vertex: */

    template<std::size_t N,class model_t>class evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);

	    /* Group type definition: */

	    DEFINE_GROUP_TYPE(model_t,SU<N>);
	    
	    /* Helper class type definition: */
	    
	    typedef colour_tensor::ff_contr_helper<SU<N>,vvvv<model_t>,0,1,2,3> helper_class;

	    /* Useful constant definitions: */

	    static const std::size_t vector_size=model_t::dimension;
	    static const std::size_t gluon_size=vector_size*(N*N-1);

	public:

            /* Initialisation phase: */

	    static void initialise()
	    {
		helper_class::initialise();
		spacetime_type::initialise();
	    }

	    /* Checking function returning the index ranges of interacting
	     * subamplitudes: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n<4)
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
		CAMGEN_ERROR_IF((A_3.range()<(N*N-1)*vector_size),"iterator 3 out of range");

		iterator B_0=A_0;
		iterator B_1=A_1;
		iterator B_2=A_2;
		iterator B_3=A_3;
		value_type C;
		
		/* Storing all subtensor inner products in memory: */
		
		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<a;++b)
		    {
			c1[a][b]=spacetime_type::dot(A_1,B_2)+spacetime_type::dot(B_1,A_2);
			c2[a][b]=spacetime_type::dot(A_1,B_3)+spacetime_type::dot(B_1,A_3);
			c3[a][b]=spacetime_type::dot(A_2,B_3)+spacetime_type::dot(B_2,A_3);
			B_1+=vector_size;
			B_2+=vector_size;
			B_3+=vector_size;
		    }
		    c1[a][a]=spacetime_type::dot(A_1,A_2);
		    c2[a][a]=spacetime_type::dot(A_1,A_3);
		    c3[a][a]=spacetime_type::dot(A_2,A_3);

		    A_1+=vector_size;
		    A_2+=vector_size;
		    A_3+=vector_size;
		    B_1-=(a*vector_size);
		    B_2-=(a*vector_size);
		    B_3-=(a*vector_size);
		}
		A_1-=gluon_size;
		A_2-=gluon_size;
		A_3-=gluon_size;

		/* Contracting with the colour tensor: */

		for(size_type a=0;a<helper_class::size();++a)
		{
		    A_0+=helper_class::jump(0,a);
		    A_1+=helper_class::jump(0,a);
		    A_2+=helper_class::jump(0,a);
		    A_3+=helper_class::jump(0,a);
		    B_0+=helper_class::jump(1,a);
		    B_1+=helper_class::jump(1,a);
		    B_2+=helper_class::jump(1,a);
		    B_3+=helper_class::jump(1,a);
		    size_type c=helper_class::index(2,a);
		    size_type d=helper_class::index(3,a);
		    C=factor*C_0*helper_class::value(a);
		    
		    /* Spacetime vertex part: */
		    
		    for(size_type mu=0;mu<vector_size;++mu)
		    {
			A_0[mu]+=C*(c2[c][d]*B_2[mu]+c1[c][d]*B_3[mu]+c3[c][d]*B_1[mu]);
			if(helper_class::index(0,a)!=helper_class::index(1,a))
			{
			    B_0[mu]+=C*(c2[c][d]*A_2[mu]+c1[c][d]*A_3[mu]+c3[c][d]*A_1[mu]);
			}
		    }
		}
		A_0-=helper_class::reset(0);
		A_1-=helper_class::reset(0);
		A_2-=helper_class::reset(0);
		A_3-=helper_class::reset(0);

	    }
	    static void second(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<(N*N-1)*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<(N*N-1)*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<(N*N-1)*vector_size),"iterator 2 out of range");
		CAMGEN_ERROR_IF((A_3.range()<(N*N-1)*vector_size),"iterator 3 out of range");
		
		iterator B_0=A_0;
		iterator B_1=A_1;
		iterator B_2=A_2;
		iterator B_3=A_3;
		value_type C;
		
		/* Storing all subtensor inner products in memory: */
		
		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<a;++b)
		    {
			c1[a][b]=spacetime_type::dot(A_0,B_2)+spacetime_type::dot(B_0,A_2);
			c2[a][b]=spacetime_type::dot(A_0,B_3)+spacetime_type::dot(B_0,A_3);
			c3[a][b]=spacetime_type::dot(A_2,B_3)+spacetime_type::dot(B_2,A_3);
			B_0+=vector_size;
			B_2+=vector_size;
			B_3+=vector_size;
		    }
		    c1[a][a]=spacetime_type::dot(A_0,A_2);
		    c2[a][a]=spacetime_type::dot(A_0,A_3);
		    c3[a][a]=spacetime_type::dot(A_2,A_3);

		    A_0+=vector_size;
		    A_2+=vector_size;
		    A_3+=vector_size;
		    B_0-=(a*vector_size);
		    B_2-=(a*vector_size);
		    B_3-=(a*vector_size);
		}
		A_0-=gluon_size;
		A_2-=gluon_size;
		A_3-=gluon_size;

		/* Contracting with the colour tensor: */

		for(size_type a=0;a<helper_class::size();++a)
		{
		    A_0+=helper_class::jump(0,a);
		    A_1+=helper_class::jump(0,a);
		    A_2+=helper_class::jump(0,a);
		    A_3+=helper_class::jump(0,a);
		    B_0+=helper_class::jump(1,a);
		    B_1+=helper_class::jump(1,a);
		    B_2+=helper_class::jump(1,a);
		    B_3+=helper_class::jump(1,a);
		    size_type c=helper_class::index(2,a);
		    size_type d=helper_class::index(3,a);
		    C=factor*C_0*helper_class::value(a);
		    
		    /* Spacetime vertex part: */
		    
		    for(size_type mu=0;mu<vector_size;++mu)
		    {
			A_1[mu]+=C*(c2[c][d]*B_2[mu]+c1[c][d]*B_3[mu]+c3[c][d]*B_0[mu]);
			if(helper_class::index(0,a)!=helper_class::index(1,a))
			{
			    B_1[mu]+=C*(c2[c][d]*A_2[mu]+c1[c][d]*A_3[mu]+c3[c][d]*A_0[mu]);
			}
		    }
		}
		A_0-=helper_class::reset(0);
		A_1-=helper_class::reset(0);
		A_2-=helper_class::reset(0);
		A_3-=helper_class::reset(0);

	    }
	    static void third(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<(N*N-1)*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<(N*N-1)*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<(N*N-1)*vector_size),"iterator 2 out of range");
		CAMGEN_ERROR_IF((A_3.range()<(N*N-1)*vector_size),"iterator 3 out of range");
		
		iterator B_0=A_0;
		iterator B_1=A_1;
		iterator B_2=A_2;
		iterator B_3=A_3;
		value_type C;
		
		/* Storing all subtensor inner products in memory: */
		
		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<a;++b)
		    {
			c1[a][b]=spacetime_type::dot(A_0,B_1)+spacetime_type::dot(B_0,A_1);
			c2[a][b]=spacetime_type::dot(A_0,B_3)+spacetime_type::dot(B_0,A_3);
			c3[a][b]=spacetime_type::dot(A_1,B_3)+spacetime_type::dot(B_1,A_3);
			B_0+=vector_size;
			B_1+=vector_size;
			B_3+=vector_size;
		    }
		    c1[a][a]=spacetime_type::dot(A_0,A_1);
		    c2[a][a]=spacetime_type::dot(A_0,A_3);
		    c3[a][a]=spacetime_type::dot(A_1,A_3);

		    A_0+=vector_size;
		    A_1+=vector_size;
		    A_3+=vector_size;
		    B_0-=(a*vector_size);
		    B_1-=(a*vector_size);
		    B_3-=(a*vector_size);
		}
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_3-=gluon_size;

		/* Contracting with the colour tensor: */

		for(size_type a=0;a<helper_class::size();++a)
		{
		    A_0+=helper_class::jump(0,a);
		    A_1+=helper_class::jump(0,a);
		    A_2+=helper_class::jump(0,a);
		    A_3+=helper_class::jump(0,a);
		    B_0+=helper_class::jump(1,a);
		    B_1+=helper_class::jump(1,a);
		    B_2+=helper_class::jump(1,a);
		    B_3+=helper_class::jump(1,a);
		    size_type c=helper_class::index(2,a);
		    size_type d=helper_class::index(3,a);
		    C=factor*C_0*helper_class::value(a);
		    
		    /* Spacetime vertex part: */
		    
		    for(size_type mu=0;mu<vector_size;++mu)
		    {
			A_2[mu]+=C*(c2[c][d]*B_1[mu]+c1[c][d]*B_3[mu]+c3[c][d]*B_0[mu]);
			if(helper_class::index(0,a)!=helper_class::index(1,a))
			{
			    B_2[mu]+=C*(c2[c][d]*A_1[mu]+c1[c][d]*A_3[mu]+c3[c][d]*A_0[mu]);
			}
		    }
		}
		A_0-=helper_class::reset(0);
		A_1-=helper_class::reset(0);
		A_2-=helper_class::reset(0);
		A_3-=helper_class::reset(0);

	    }
	    static void fourth(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<(N*N-1)*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<(N*N-1)*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<(N*N-1)*vector_size),"iterator 2 out of range");
		CAMGEN_ERROR_IF((A_3.range()<(N*N-1)*vector_size),"iterator 3 out of range");
		
		iterator B_0=A_0;
		iterator B_1=A_1;
		iterator B_2=A_2;
		iterator B_3=A_3;
		value_type C;
		
		/* Storing all subtensor inner products in memory: */
		
		for(size_type a=0;a<N*N-1;++a)
		{
		    for(size_type b=0;b<a;++b)
		    {
			c1[a][b]=spacetime_type::dot(A_0,B_1)+spacetime_type::dot(B_0,A_1);
			c2[a][b]=spacetime_type::dot(A_0,B_2)+spacetime_type::dot(B_0,A_2);
			c3[a][b]=spacetime_type::dot(A_1,B_2)+spacetime_type::dot(B_1,A_2);
			B_0+=vector_size;
			B_1+=vector_size;
			B_2+=vector_size;
		    }
		    c1[a][a]=spacetime_type::dot(A_0,A_1);
		    c2[a][a]=spacetime_type::dot(A_0,A_2);
		    c3[a][a]=spacetime_type::dot(A_1,A_2);

		    A_0+=vector_size;
		    A_1+=vector_size;
		    A_2+=vector_size;
		    B_0-=(a*vector_size);
		    B_1-=(a*vector_size);
		    B_2-=(a*vector_size);
		}
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_2-=gluon_size;

		/* Contracting with the colour tensor: */

		for(size_type a=0;a<helper_class::size();++a)
		{
		    A_0+=helper_class::jump(0,a);
		    A_1+=helper_class::jump(0,a);
		    A_2+=helper_class::jump(0,a);
		    A_3+=helper_class::jump(0,a);
		    B_0+=helper_class::jump(1,a);
		    B_1+=helper_class::jump(1,a);
		    B_2+=helper_class::jump(1,a);
		    B_3+=helper_class::jump(1,a);
		    size_type c=helper_class::index(2,a);
		    size_type d=helper_class::index(3,a);
		    C=factor*C_0*helper_class::value(a);
		    
		    /* Spacetime vertex part: */
		    
		    for(size_type mu=0;mu<vector_size;++mu)
		    {
			A_3[mu]+=C*(c2[c][d]*B_1[mu]+c1[c][d]*B_2[mu]+c3[c][d]*B_0[mu]);
			if(helper_class::index(0,a)!=helper_class::index(1,a))
			{
			    B_3[mu]+=C*(c2[c][d]*A_1[mu]+c1[c][d]*A_2[mu]+c3[c][d]*A_0[mu]);
			}
		    }
		}
		A_0-=helper_class::reset(0);
		A_1-=helper_class::reset(0);
		A_2-=helper_class::reset(0);
		A_3-=helper_class::reset(0);

	    }
	private:

	    /* Inner product data holders: */

	    static value_type c1[N*N-1][N*N-1];
	    static value_type c2[N*N-1][N*N-1];
	    static value_type c3[N*N-1][N*N-1];

	    /* Initialisation tag: */

	    static bool initialised;
    };	
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::vector_size;
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::gluon_size;
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::c1[N*N-1][N*N-1]={{0}};
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::c2[N*N-1][N*N-1]={{0}};
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::c3[N*N-1][N*N-1]={{0}};
    template<std::size_t N,class model_t>bool evaluate<compose_vertices<colour_tensor::ff_contr<SU<N>,0,1,2,3>,vvvv<model_t> > >::initialised=false;

    /* Specialisation of the 4-gluon vertex with the gluons in the colour-flow
     * representation: */

    template<std::size_t N,class model_t>class compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >
    {
	private:

	    /* Useful constant definition: */

	    static const std::size_t gluon_size=N*N*(model_t::dimension);

	public:
	    
	    /* Refernce definition to the model type: */
	    
	    typedef model_t model_type;

	    /* Rank of the vertex: */

	    static const std::size_t rank=4;
	    
	    /* Number of couplings used by the vertex: */
	    
	    static const std::size_t params=1;

	    /* Boolean denoting whether the Feynman rule depends on momenta: */

	    static const bool p_dependent=false;
	    
	    /* Boolean denoting whether the Feynman rule is fermionic: */
	    
	    static const bool fermionic=false;

	    /* Size of the vertex tensor: */

	    static const std::size_t tensor_size=gluon_size*gluon_size*gluon_size*gluon_size;

	    /* Subamplitude sizes: */

	    static const std::size_t sizes[4];

	    /* Function returning the formula of the Feynman rule: */

	    static const std::string formula;
    };
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::gluon_size;
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::rank;
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::params;
    template<std::size_t N,class model_t>const bool compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::p_dependent;
    template<std::size_t N,class model_t>const bool compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::fermionic;
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::tensor_size;
    template<std::size_t N,class model_t>const std::size_t compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::sizes[4]={N*N*model_t::dimension,N*N*model_t::dimension,N*N*model_t::dimension,N*N*model_t::dimension};
    template<std::size_t N,class model_t>const std::string compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> >::formula="1/8d(A1,B2)d(A2,B3)d(A3,B4)d(A4,B1)(2g(mu1,mu3)g(mu2,mu4)-g(mu1,mu2)g(mu3,mu4)-g(mu1,mu4)g(mu2,mu3))+2<->3<->4";

    /* Specialisation of the evaluate class template for the gluon 4-vertex with
     * gluons in the colour-flow representation: */

    template<std::size_t N,class model_t>class evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > vertex_type;
	    
	    /* The usual type definitions: */
	    
	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);
	    
	    /* Group type definition: */
	    
	    DEFINE_GROUP_TYPE(model_t,SU<N>);

	    /* Some useful constants: */

	    static const std::size_t vector_size=model_t::dimension;
	    static const std::size_t quark_size=N*vector_size;
	    static const std::size_t gluon_size=N*quark_size;

	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
	    }

	    /* Checking function returning the index ranges of the interacting
	     * gluons: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n<4)
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
		CAMGEN_ERROR_IF((A_3.range()<N*N*vector_size),"iterator 3 out of range");

		/* Storing inner products in memory: */

		for(size_type a=0;a<N;++a)
		{
		    for(size_type b=0;b<N;++b)
		    {
			c12[a][b]=0;
			c13[a][b]=0;
			c23[a][b]=0;
		    }
		}
		iterator B_1=A_1;
		iterator B_2=A_2;
		iterator B_3=A_3;
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    c12[a][b]+=spacetime_type::dot(A_1,B_2);
			    c12[a][b]+=spacetime_type::dot(B_1,A_2);
			    c13[a][b]+=spacetime_type::dot(A_1,B_3);
			    c13[a][b]+=spacetime_type::dot(B_1,A_3);
			    c23[a][b]+=spacetime_type::dot(A_2,B_3);
			    c23[a][b]+=spacetime_type::dot(B_2,A_3);

			    B_1+=vector_size;
			    B_2+=vector_size;
			    B_3+=vector_size;
			}
			A_1+=vector_size;
			A_2+=vector_size;
			A_3+=vector_size;
		    }
		    B_1-=gluon_size;
		    B_2-=gluon_size;
		    B_3-=gluon_size;
		}
		A_1-=gluon_size;
		A_2-=gluon_size;
		A_3-=gluon_size;
		iterator D_1=A_1;
		iterator D_2=A_2;
		iterator D_3=A_3;
		iterator E_1=A_1;
		iterator E_2=A_2;
		iterator E_3=A_3;
		value_type temp1;
		value_type temp2;
		value_type temp3;

		/* In the inner loop: A_0 is G_0(ab), A_i is G_i(ac), B_i is
		 * G_i(cb), D_i is G_i(cd), and E_i is G_i(db) */
		
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    for(size_type mu=0;mu<vector_size;++mu)
			    {
				A_0[mu]+=factor*C_0*(A_1[mu]*c23[c][b]+c23[a][c]*B_1[mu]+A_2[mu]*c13[c][b]+c13[a][c]*B_2[mu]+A_3[mu]*c12[c][b]+c12[a][c]*B_3[mu]);
			    }	
			    for(size_type d=0;d<N;++d)
			    {
				temp1=spacetime_type::dot(A_2,E_3)+spacetime_type::dot(A_3,E_2);
				temp2=spacetime_type::dot(A_1,E_3)+spacetime_type::dot(A_3,E_1);
				temp3=spacetime_type::dot(A_1,E_2)+spacetime_type::dot(A_2,E_1);

				for(size_type mu=0;mu<vector_size;++mu)
				{
				    A_0[mu]-=(r_value_type)2*factor*C_0*(D_1[mu]*temp1+D_2[mu]*temp2+D_3[mu]*temp3);
				}
				D_1+=vector_size;
				D_2+=vector_size;
				D_3+=vector_size;
				E_1+=quark_size;
				E_2+=quark_size;
				E_3+=quark_size;
			    }			
			    A_0+=vector_size;
			    B_1+=vector_size;
			    B_2+=vector_size;
			    B_3+=vector_size;
			    D_1-=quark_size;
			    D_2-=quark_size;
			    D_3-=quark_size;
			    E_1-=(gluon_size-vector_size);
			    E_2-=(gluon_size-vector_size);
			    E_3-=(gluon_size-vector_size);
			}
			A_0-=quark_size;
			A_1+=vector_size;
			A_2+=vector_size;
			A_3+=vector_size;
			D_1+=quark_size;
			D_2+=quark_size;
			D_3+=quark_size;
			E_1-=quark_size;
			E_2-=quark_size;
			E_3-=quark_size;
		    }
		    A_0+=quark_size;
		    B_1-=gluon_size;
		    B_2-=gluon_size;
		    B_3-=gluon_size;
		    D_1-=gluon_size;
		    D_2-=gluon_size;
		    D_3-=gluon_size;
		}
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_2-=gluon_size;
		A_3-=gluon_size;

	    }
	    static void second(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<N*N*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<N*N*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<N*N*vector_size),"iterator 2 out of range");
		CAMGEN_ERROR_IF((A_3.range()<N*N*vector_size),"iterator 3 out of range");

		/* Storing inner products in memory: */

		for(size_type a=0;a<N;++a)
		{
		    for(size_type b=0;b<N;++b)
		    {
			c12[a][b]=0;
			c13[a][b]=0;
			c23[a][b]=0;
		    }
		}
		iterator B_0=A_0;
		iterator B_2=A_2;
		iterator B_3=A_3;
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    c12[a][b]+=spacetime_type::dot(A_0,B_2);
			    c12[a][b]+=spacetime_type::dot(B_0,A_2);
			    c13[a][b]+=spacetime_type::dot(A_0,B_3);
			    c13[a][b]+=spacetime_type::dot(B_0,A_3);
			    c23[a][b]+=spacetime_type::dot(A_2,B_3);
			    c23[a][b]+=spacetime_type::dot(B_2,A_3);

			    B_0+=vector_size;
			    B_2+=vector_size;
			    B_3+=vector_size;
			}
			A_0+=vector_size;
			A_2+=vector_size;
			A_3+=vector_size;
		    }
		    B_0-=gluon_size;
		    B_2-=gluon_size;
		    B_3-=gluon_size;
		}
		A_0-=gluon_size;
		A_2-=gluon_size;
		A_3-=gluon_size;
		iterator D_0=A_0;
		iterator D_2=A_2;
		iterator D_3=A_3;
		iterator E_0=A_0;
		iterator E_2=A_2;
		iterator E_3=A_3;
		value_type temp1;
		value_type temp2;
		value_type temp3;
		
		/* In the inner loop: A_1 is G_0(ab), A_i is G_i(ac), B_i is
		 * G_i(cb), D_i is G_i(cd), and E_i is G_i(db): */
		
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    for(size_type mu=0;mu<vector_size;++mu)
			    {
				A_1[mu]+=factor*C_0*(A_0[mu]*c23[c][b]+c23[a][c]*B_0[mu]+A_2[mu]*c13[c][b]+c13[a][c]*B_2[mu]+A_3[mu]*c12[c][b]+c12[a][c]*B_3[mu]);
			    }	
			    for(size_type d=0;d<N;++d)
			    {
				temp1=spacetime_type::dot(A_2,E_3)+spacetime_type::dot(A_3,E_2);
				temp2=spacetime_type::dot(A_0,E_3)+spacetime_type::dot(A_3,E_0);
				temp3=spacetime_type::dot(A_0,E_2)+spacetime_type::dot(A_2,E_0);

				for(size_type mu=0;mu<vector_size;++mu)
				{
				    A_1[mu]-=(r_value_type)2*factor*C_0*(D_0[mu]*temp1+D_2[mu]*temp2+D_3[mu]*temp3);
				}
				D_0+=vector_size;
				D_2+=vector_size;
				D_3+=vector_size;
				E_0+=quark_size;
				E_2+=quark_size;
				E_3+=quark_size;
			    }			
			    A_1+=vector_size;
			    B_0+=vector_size;
			    B_2+=vector_size;
			    B_3+=vector_size;
			    D_0-=quark_size;
			    D_2-=quark_size;
			    D_3-=quark_size;
			    E_0-=(gluon_size-vector_size);
			    E_2-=(gluon_size-vector_size);
			    E_3-=(gluon_size-vector_size);
			}
			A_1-=quark_size;
			A_0+=vector_size;
			A_2+=vector_size;
			A_3+=vector_size;
			D_0+=quark_size;
			D_2+=quark_size;
			D_3+=quark_size;
			E_0-=quark_size;
			E_2-=quark_size;
			E_3-=quark_size;
		    }
		    A_1+=quark_size;
		    B_0-=gluon_size;
		    B_2-=gluon_size;
		    B_3-=gluon_size;
		    D_0-=gluon_size;
		    D_2-=gluon_size;
		    D_3-=gluon_size;
		}
		A_1-=gluon_size;
		A_0-=gluon_size;
		A_2-=gluon_size;
		A_3-=gluon_size;

	    }
	    static void third(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<N*N*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<N*N*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<N*N*vector_size),"iterator 2 out of range");
		CAMGEN_ERROR_IF((A_3.range()<N*N*vector_size),"iterator 3 out of range");

		/* Storing inner products in memory: */

		for(size_type a=0;a<N;++a)
		{
		    for(size_type b=0;b<N;++b)
		    {
			c12[a][b]=0;
			c13[a][b]=0;
			c23[a][b]=0;
		    }
		}
		iterator B_0=A_0;
		iterator B_1=A_1;
		iterator B_3=A_3;
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    c12[a][b]+=spacetime_type::dot(A_0,B_1);
			    c12[a][b]+=spacetime_type::dot(B_0,A_1);
			    c13[a][b]+=spacetime_type::dot(A_0,B_3);
			    c13[a][b]+=spacetime_type::dot(B_0,A_3);
			    c23[a][b]+=spacetime_type::dot(A_1,B_3);
			    c23[a][b]+=spacetime_type::dot(B_1,A_3);

			    B_0+=vector_size;
			    B_1+=vector_size;
			    B_3+=vector_size;
			}
			A_0+=vector_size;
			A_1+=vector_size;
			A_3+=vector_size;
		    }
		    B_0-=gluon_size;
		    B_1-=gluon_size;
		    B_3-=gluon_size;
		}
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_3-=gluon_size;
		iterator D_0=A_0;
		iterator D_1=A_1;
		iterator D_3=A_3;
		iterator E_0=A_0;
		iterator E_1=A_1;
		iterator E_3=A_3;
		value_type temp1;
		value_type temp2;
		value_type temp3;
		
		/* In the inner loop: A_2 is G_0(ab), A_i is G_i(ac), B_i is
		 * G_i(cb), D_i is G_i(cd), and E_i is G_i(db) */
		
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    for(size_type mu=0;mu<vector_size;++mu)
			    {
				A_2[mu]+=factor*C_0*(A_0[mu]*c23[c][b]+c23[a][c]*B_0[mu]+A_1[mu]*c13[c][b]+c13[a][c]*B_1[mu]+A_3[mu]*c12[c][b]+c12[a][c]*B_3[mu]);
			    }	
			    for(size_type d=0;d<N;++d)
			    {
				temp1=spacetime_type::dot(A_1,E_3)+spacetime_type::dot(A_3,E_1);
				temp2=spacetime_type::dot(A_0,E_3)+spacetime_type::dot(A_3,E_0);
				temp3=spacetime_type::dot(A_0,E_1)+spacetime_type::dot(A_1,E_0);

				for(size_type mu=0;mu<vector_size;++mu)
				{
				    A_2[mu]-=(r_value_type)2*factor*C_0*(D_0[mu]*temp1+D_1[mu]*temp2+D_3[mu]*temp3);
				}
				D_0+=vector_size;
				D_1+=vector_size;
				D_3+=vector_size;
				E_0+=quark_size;
				E_1+=quark_size;
				E_3+=quark_size;
			    }
			    A_2+=vector_size;
			    B_0+=vector_size;
			    B_1+=vector_size;
			    B_3+=vector_size;
			    D_0-=quark_size;
			    D_1-=quark_size;
			    D_3-=quark_size;
			    E_0-=(gluon_size-vector_size);
			    E_1-=(gluon_size-vector_size);
			    E_3-=(gluon_size-vector_size);
			}
			A_2-=quark_size;
			A_0+=vector_size;
			A_1+=vector_size;
			A_3+=vector_size;
			D_0+=quark_size;
			D_1+=quark_size;
			D_3+=quark_size;
			E_0-=quark_size;
			E_1-=quark_size;
			E_3-=quark_size;
		    }
		    A_2+=quark_size;
		    B_0-=gluon_size;
		    B_1-=gluon_size;
		    B_3-=gluon_size;
		    D_0-=gluon_size;
		    D_1-=gluon_size;
		    D_3-=gluon_size;
		}
		A_2-=gluon_size;
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_3-=gluon_size;

	    }
	    static void fourth(ARG_LIST)
	    {
		
		CAMGEN_ERROR_IF((A_0.range()<N*N*vector_size),"iterator 0 out of range");
		CAMGEN_ERROR_IF((A_1.range()<N*N*vector_size),"iterator 1 out of range");
		CAMGEN_ERROR_IF((A_2.range()<N*N*vector_size),"iterator 2 out of range");
		CAMGEN_ERROR_IF((A_3.range()<N*N*vector_size),"iterator 3 out of range");

		/* Storing inner products in memory: */

		for(size_type a=0;a<N;++a)
		{
		    for(size_type b=0;b<N;++b)
		    {
			c12[a][b]=0;
			c13[a][b]=0;
			c23[a][b]=0;
		    }
		}
		iterator B_0=A_0;
		iterator B_1=A_1;
		iterator B_2=A_2;
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    c12[a][b]+=spacetime_type::dot(A_0,B_1);
			    c12[a][b]+=spacetime_type::dot(B_0,A_1);
			    c13[a][b]+=spacetime_type::dot(A_0,B_2);
			    c13[a][b]+=spacetime_type::dot(B_0,A_2);
			    c23[a][b]+=spacetime_type::dot(A_1,B_2);
			    c23[a][b]+=spacetime_type::dot(B_1,A_2);

			    B_0+=vector_size;
			    B_1+=vector_size;
			    B_2+=vector_size;
			}
			A_0+=vector_size;
			A_1+=vector_size;
			A_2+=vector_size;
		    }
		    B_0-=gluon_size;
		    B_1-=gluon_size;
		    B_2-=gluon_size;
		}
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_2-=gluon_size;
		iterator D_0=A_0;
		iterator D_1=A_1;
		iterator D_2=A_2;
		iterator E_0=A_0;
		iterator E_1=A_1;
		iterator E_2=A_2;
		value_type temp1;
		value_type temp2;
		value_type temp3;

		/* In the inner loop: A_3 is G_0(ab), A_i is G_i(ac), B_i is
		 * G_i(cb), D_i is G_i(cd), and E_i is G_i(db) */
		
		for(size_type a=0;a<N;++a)
		{
		    for(size_type c=0;c<N;++c)
		    {
			for(size_type b=0;b<N;++b)
			{
			    for(size_type mu=0;mu<vector_size;++mu)
			    {
				A_3[mu]+=factor*C_0*(A_0[mu]*c23[c][b]+c23[a][c]*B_0[mu]+A_1[mu]*c13[c][b]+c13[a][c]*B_1[mu]+A_2[mu]*c12[c][b]+c12[a][c]*B_2[mu]);
			    }	
			    for(size_type d=0;d<N;++d)
			    {
				temp1=spacetime_type::dot(A_1,E_2)+spacetime_type::dot(A_2,E_1);
				temp2=spacetime_type::dot(A_0,E_2)+spacetime_type::dot(A_2,E_0);
				temp3=spacetime_type::dot(A_0,E_1)+spacetime_type::dot(A_1,E_0);

				for(size_type mu=0;mu<vector_size;++mu)
				{
				    A_3[mu]-=(r_value_type)2*factor*C_0*(D_0[mu]*temp1+D_1[mu]*temp2+D_2[mu]*temp3);
				}
				D_0+=vector_size;
				D_1+=vector_size;
				D_2+=vector_size;
				E_0+=quark_size;
				E_1+=quark_size;
				E_2+=quark_size;
			    }
			    A_3+=vector_size;
			    B_0+=vector_size;
			    B_1+=vector_size;
			    B_2+=vector_size;
			    D_0-=quark_size;
			    D_1-=quark_size;
			    D_2-=quark_size;
			    E_0-=(gluon_size-vector_size);
			    E_1-=(gluon_size-vector_size);
			    E_2-=(gluon_size-vector_size);
			}
			A_3-=quark_size;
			A_0+=vector_size;
			A_1+=vector_size;
			A_2+=vector_size;
			D_0+=quark_size;
			D_1+=quark_size;
			D_2+=quark_size;
			E_0-=quark_size;
			E_1-=quark_size;
			E_2-=quark_size;
		    }
		    A_3+=quark_size;
		    B_0-=gluon_size;
		    B_1-=gluon_size;
		    B_2-=gluon_size;
		    D_0-=gluon_size;
		    D_1-=gluon_size;
		    D_2-=gluon_size;
		}
		A_3-=gluon_size;
		A_0-=gluon_size;
		A_1-=gluon_size;
		A_2-=gluon_size;

	    }

	private:

	    /* Inner product data holders: */

	    static value_type c12[N][N];
	    static value_type c13[N][N];
	    static value_type c23[N][N];

	    /* Initialisation tag: */

	    static bool initialised;
    };
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::vector_size;
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::quark_size;
    template<std::size_t N,class model_t>const std::size_t evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::gluon_size;
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::c12[N][N]={{0}};
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::c13[N][N]={{0}};
    template<std::size_t N,class model_t>typename evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::value_type evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::c23[N][N]={{0}};
    template<std::size_t N,class model_t>bool evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::initialised=false;

    /* Specialisation of the cfd_evaluate class template for the gluon 4-vertex
     * with gluons in the colour-flow representation: */

    template<std::size_t N,class model_t>class cfd_evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >
    {
	public:

	    /* Vertex type definition: */

	    typedef compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > vertex_type;

	    /* The usual type definitions: */

	    DEFINE_BASIC_TYPES(model_t);

	private:

	    /* Spacetime type definition: */

	    DEFINE_SPACETIME_TYPE(model_t);
	    
	    /* Group type definition: */
	    
	    DEFINE_GROUP_TYPE(model_t,SU<N>);

	    /* Some useful constants: */

	    static const std::size_t vector_size=model_t::dimension;
	
	public:

	    /* Initialisation phase: */

	    static void initialise()
	    {
		spacetime_type::initialise();
	    }
	    
	    /* Checking function returning the index ranges of the interacting
	     * gluons: */

	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		v.clear();
		if(n<4)
		{
		    v.push_back(N);
		    v.push_back(N);
		    v.push_back(vector_size);
		}
		return v;
	    }
	    
	    /* Recursive relations: */

	    static void first(CFD_ARG_LIST)
	    {

		size_type S1=(A_1.get_offset()/vector_size)%(N*N);
		size_type A1=S1/N;
		size_type B1=S1%N;
		size_type S2=(A_2.get_offset()/vector_size)%(N*N);
		size_type A2=S2/N;
		size_type B2=S2%N;
		size_type S3=(A_3.get_offset()/vector_size)%(N*N);
		size_type A3=S3/N;
		size_type B3=S3%N;

		if(B1==A2 and B2==A3)
		{
		    A_0+=(A1*N+B3)*vector_size;
		    apply(factor*C_0,A_0,A_1,A_2,A_3);
		    produced_iters.insert(A_0);
		    A_0-=(A1*N+B3)*vector_size;
		}
		if(B1==A3 and B3==A2)
		{
		    A_0+=(A1*N+B2)*vector_size;
		    apply(factor*C_0,A_0,A_1,A_3,A_2);
		    produced_iters.insert(A_0);
		    A_0-=(A1*N+B2)*vector_size;
		}
		if(B2==A1 and B1==A3)
		{
		    A_0+=(A2*N+B3)*vector_size;
		    apply(factor*C_0,A_0,A_2,A_1,A_3);
		    produced_iters.insert(A_0);
		    A_0-=(A2*N+B3)*vector_size;
		}
		if(B2==A3 and B3==A1)
		{
		    A_0+=(A2*N+B1)*vector_size;
		    apply(factor*C_0,A_0,A_2,A_3,A_1);
		    produced_iters.insert(A_0);
		    A_0-=(A2*N+B1)*vector_size;
		}
		if(B3==A1 and B1==A2)
		{
		    A_0+=(A3*N+B2)*vector_size;
		    apply(factor*C_0,A_0,A_3,A_1,A_2);
		    produced_iters.insert(A_0);
		    A_0-=(A3*N+B2)*vector_size;
		}
		if(B3==A2 and B2==A1)
		{
		    A_0+=(A3*N+B1)*vector_size;
		    apply(factor*C_0,A_0,A_3,A_2,A_1);
		    produced_iters.insert(A_0);
		    A_0-=(A3*N+B1)*vector_size;
		}
	    }
	    static void second(CFD_ARG_LIST)
	    {

		size_type S0=(A_0.get_offset()/vector_size)%(N*N);
		size_type A0=S0/N;
		size_type B0=S0%N;
		size_type S2=(A_2.get_offset()/vector_size)%(N*N);
		size_type A2=S2/N;
		size_type B2=S2%N;
		size_type S3=(A_3.get_offset()/vector_size)%(N*N);
		size_type A3=S3/N;
		size_type B3=S3%N;

		if(B0==A2 and B2==A3)
		{
		    A_1+=(A0*N+B3)*vector_size;
		    apply(factor*C_0,A_1,A_0,A_2,A_3);
		    produced_iters.insert(A_1);
		    A_1-=(A0*N+B3)*vector_size;
		}
		if(B0==A3 and B3==A2)
		{
		    A_1+=(A0*N+B2)*vector_size;
		    apply(factor*C_0,A_1,A_0,A_3,A_2);
		    produced_iters.insert(A_1);
		    A_1-=(A0*N+B2)*vector_size;
		}
		if(B2==A0 and B0==A3)
		{
		    A_1+=(A2*N+B3)*vector_size;
		    apply(factor*C_0,A_1,A_2,A_0,A_3);
		    produced_iters.insert(A_1);
		    A_1-=(A2*N+B3)*vector_size;
		}
		if(B2==A3 and B3==A0)
		{
		    A_1+=(A2*N+B0)*vector_size;
		    apply(factor*C_0,A_1,A_2,A_3,A_0);
		    produced_iters.insert(A_1);
		    A_1-=(A2*N+B0)*vector_size;
		}
		if(B3==A0 and B0==A2)
		{
		    A_1+=(A3*N+B2)*vector_size;
		    apply(factor*C_0,A_1,A_3,A_0,A_2);
		    produced_iters.insert(A_1);
		    A_1-=(A3*N+B2)*vector_size;
		}
		if(B3==A2 and B2==A0)
		{
		    A_1+=(A3*N+B0)*vector_size;
		    apply(factor*C_0,A_1,A_3,A_2,A_0);
		    produced_iters.insert(A_1);
		    A_1-=(A3*N+B0)*vector_size;
		}
	    }
	    static void third(CFD_ARG_LIST)
	    {

		size_type S0=(A_0.get_offset()/vector_size)%(N*N);
		size_type A0=S0/N;
		size_type B0=S0%N;
		size_type S1=(A_1.get_offset()/vector_size)%(N*N);
		size_type A1=S1/N;
		size_type B1=S1%N;
		size_type S3=(A_3.get_offset()/vector_size)%(N*N);
		size_type A3=S3/N;
		size_type B3=S3%N;

		if(B1==A0 and B0==A3)
		{
		    A_2+=(A1*N+B3)*vector_size;
		    apply(factor*C_0,A_2,A_1,A_0,A_3);
		    produced_iters.insert(A_2);
		    A_2-=(A1*N+B3)*vector_size;
		}
		if(B1==A3 and B3==A0)
		{
		    A_2+=(A1*N+B0)*vector_size;
		    apply(factor*C_0,A_2,A_1,A_3,A_0);
		    produced_iters.insert(A_2);
		    A_2-=(A1*N+B0)*vector_size;
		}
		if(B0==A1 and B1==A3)
		{
		    A_2+=(A0*N+B3)*vector_size;
		    apply(factor*C_0,A_2,A_0,A_1,A_3);
		    produced_iters.insert(A_2);
		    A_2-=(A0*N+B3)*vector_size;
		}
		if(B0==A3 and B3==A1)
		{
		    A_2+=(A0*N+B1)*vector_size;
		    apply(factor*C_0,A_2,A_0,A_3,A_1);
		    produced_iters.insert(A_2);
		    A_2-=(A0*N+B1)*vector_size;
		}
		if(B3==A1 and B1==A0)
		{
		    A_2+=(A3*N+B0)*vector_size;
		    apply(factor*C_0,A_2,A_3,A_1,A_0);
		    produced_iters.insert(A_2);
		    A_2-=(A3*N+B0)*vector_size;
		}
		if(B3==A0 and B0==A1)
		{
		    A_2+=(A3*N+B1)*vector_size;
		    apply(factor*C_0,A_2,A_3,A_0,A_1);
		    produced_iters.insert(A_2);
		    A_2-=(A3*N+B1)*vector_size;
		}
	    }
	    static void fourth(CFD_ARG_LIST)
	    {

		size_type S0=(A_0.get_offset()/vector_size)%(N*N);
		size_type A0=S0/N;
		size_type B0=S0%N;
		size_type S1=(A_1.get_offset()/vector_size)%(N*N);
		size_type A1=S1/N;
		size_type B1=S1%N;
		size_type S2=(A_2.get_offset()/vector_size)%(N*N);
		size_type A2=S2/N;
		size_type B2=S2%N;

		if(B1==A2 and B2==A0)
		{
		    A_3+=(A1*N+B0)*vector_size;
		    apply(factor*C_0,A_3,A_1,A_2,A_0);
		    produced_iters.insert(A_3);
		    A_3-=(A1*N+B0)*vector_size;
		}
		if(B1==A0 and B0==A2)
		{
		    A_3+=(A1*N+B2)*vector_size;
		    apply(factor*C_0,A_3,A_1,A_0,A_2);
		    produced_iters.insert(A_3);
		    A_3-=(A1*N+B2)*vector_size;
		}
		if(B2==A1 and B1==A0)
		{
		    A_3+=(A2*N+B0)*vector_size;
		    apply(factor*C_0,A_3,A_2,A_1,A_0);
		    produced_iters.insert(A_3);
		    A_3-=(A2*N+B0)*vector_size;
		}
		if(B2==A0 and B0==A1)
		{
		    A_3+=(A2*N+B1)*vector_size;
		    apply(factor*C_0,A_3,A_2,A_0,A_1);
		    produced_iters.insert(A_3);
		    A_3-=(A2*N+B1)*vector_size;
		}
		if(B0==A1 and B1==A2)
		{
		    A_3+=(A0*N+B2)*vector_size;
		    apply(factor*C_0,A_3,A_0,A_1,A_2);
		    produced_iters.insert(A_3);
		    A_3-=(A0*N+B2)*vector_size;
		}
		if(B0==A2 and B2==A1)
		{
		    A_3+=(A0*N+B1)*vector_size;
		    apply(factor*C_0,A_3,A_0,A_2,A_1);
		    produced_iters.insert(A_3);
		    A_3-=(A0*N+B1)*vector_size;
		}
	    }

	private:

	    static void apply(const value_type c,iterator it0,iterator it1,iterator it2,iterator it3)
	    {
		CAMGEN_ERROR_IF((it0.range()<vector_size),"target iterator out of range");
		CAMGEN_ERROR_IF((it1.range()<vector_size),"first argument iterator out of range");
		CAMGEN_ERROR_IF((it2.range()<vector_size),"second argument iterator out of range");
		CAMGEN_ERROR_IF((it3.range()<vector_size),"third argument iterator out of range");

		value_type c1=spacetime_type::dot(it1,it2);
		value_type c2=spacetime_type::dot(it1,it3);
		value_type c3=spacetime_type::dot(it2,it3);
		for(size_type mu=0;mu<model_t::dimension;++mu)
		{
		    it0[mu]-=c*((r_value_type)2*c2*it2[mu]-c1*it3[mu]-c3*it1[mu]);
		}
	    }
    };
    template<std::size_t N,class model_t>const std::size_t cfd_evaluate<compose_vertices<colour_tensor::CF_ff_contr<N,0,1,2,3>,vvvv<model_t> > >::vector_size;
}


#endif /*CAMGEN_GGGG_H_*/
