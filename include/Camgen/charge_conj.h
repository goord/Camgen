//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_CHARGE_CONJ_H_
#define CAMGEN_CHARGE_CONJ_H_

#include <vector>
#include <Camgen/eval.h>
#include <Camgen/comp_vert.h>
#include <Camgen/comp_contr.h>
#include <Camgen/def_args.h>

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Declaration and definition of the charge conjugation class templates. The     *
 * purpose of these classes is to transform the spacetime part of fermionic      *
 * vertex and contraction classes to fermion-current reversing recursive         *
 * relations. The method is to invoke additional static 'reversed' functions in  *
 * the evaluate class, which should therefore be defined in all specialisations  *
 * of the evaluate class for fermionic spacetime structures.                     *
 *                                                                               *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{

    /* Definition of the charge conjugate contraction class, holding only
     * compile-time information identical to the original contraction class (the
     * template parameter): */

    template<class contr_t,bool ferm=contr_t::fermionic>class charge_conj_contr
    {
	public:
	    static const bool fermionic=contr_t::fermionic;
	    static const std::size_t tensor_size=contr_t::tensor_size;
	    static const std::size_t size=contr_t::size;
    };
    template<class contr_t,bool ferm>const bool charge_conj_contr<contr_t,ferm>::fermionic;
    template<class contr_t,bool ferm>const std::size_t charge_conj_contr<contr_t,ferm>::tensor_size;
    template<class contr_t,bool ferm>const std::size_t charge_conj_contr<contr_t,ferm>::size;

    /* Specialisation of the evaluate class for charge conjugated template
     * arguments. For bosonic contractions, the methods give exactly the same result 
     * as the original evaluation class, but fermion contractions will be specialised: */
    
    template<class contr_t>class evaluate< charge_conj_contr<contr_t,false> >
    {
	public:
	    typedef contr_t contraction_type;
	    typedef typename evaluate<contraction_type>::model_type model_type;
	    typedef typename evaluate<contraction_type>::value_type value_type;
	    typedef typename evaluate<contraction_type>::r_value_type r_value_type;
	    typedef typename evaluate<contraction_type>::tensor_type tensor_type;
	    typedef typename evaluate<contraction_type>::size_type size_type;
	    typedef typename evaluate<contraction_type>::iterator iterator;
	    typedef typename evaluate<contraction_type>::const_iterator const_iterator;
	    typedef typename evaluate<contraction_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<contraction_type>::initialise();
	    }
	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		return evaluate<contraction_type>::fill_rank_vector(r);
	    }
	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		return evaluate<contraction_type>::apply(it1,it2,colours,c);
	    }
	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		evaluate<contraction_type>::apply(result,it1,it2);
	    }
    };
    
    /* Specialisation of the evaluate class for charge conjugated template
     * arguments. For fermionic contractions, the methods call the Cc_apply 
     * function in the fermionic contraction evaluation: */
    
    template<class contr_t>class evaluate< charge_conj_contr<contr_t,true> >
    {
	public:
	    typedef contr_t contraction_type;
	    typedef typename evaluate<contraction_type>::model_type model_type;
	    typedef typename evaluate<contraction_type>::value_type value_type;
	    typedef typename evaluate<contraction_type>::r_value_type r_value_type;
	    typedef typename evaluate<contraction_type>::tensor_type tensor_type;
	    typedef typename evaluate<contraction_type>::size_type size_type;
	    typedef typename evaluate<contraction_type>::iterator iterator;
	    typedef typename evaluate<contraction_type>::const_iterator const_iterator;
	    typedef typename evaluate<contraction_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<contraction_type>::initialise();
	    }
	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		return evaluate<contraction_type>::fill_rank_vector(r);
	    }
	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		return evaluate<contraction_type>::Cc_apply(it1,it2,colours,c);
	    }
	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		evaluate<contraction_type>::Cc_apply(result,it1,it2);
	    }
    };

    /* Specialisation of the evaluate class for composed contractions. All
     * member functions are passed from the composition of the outer structure
     * and the charge conjugated inner contraction: */

    template<class contr_t1,class contr_t2>class evaluate<charge_conj_contr< compose_contraction<contr_t1,contr_t2>,true> >
    {
	public:
	    typedef compose_contraction<contr_t1,charge_conj_contr<contr_t2> > contraction_type;
	    typedef typename evaluate<contraction_type>::model_type model_type;
	    typedef typename evaluate<contraction_type>::value_type value_type;
	    typedef typename evaluate<contraction_type>::r_value_type r_value_type;
	    typedef typename evaluate<contraction_type>::tensor_type tensor_type;
	    typedef typename evaluate<contraction_type>::size_type size_type;
	    typedef typename evaluate<contraction_type>::iterator iterator;
	    typedef typename evaluate<contraction_type>::const_iterator const_iterator;
	    typedef typename evaluate<contraction_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<contraction_type>::initialise();
	    }
	    static std::vector<size_type>& fill_rank_vector(std::vector<size_type>& r)
	    {
		return evaluate<contraction_type>::fill_rank_vector(r);
	    }
	    static value_type apply(const_iterator it1,const_iterator it2,const std::vector<size_type>& colours,size_type c)
	    {
		return evaluate<contraction_type>::apply(it1,it2,colours,c);
	    }
	    static void apply(value_type& result,const_iterator it1,const_iterator it2)
	    {
		evaluate<contraction_type>::apply(result,it1,it2);
	    }
    };

    /* Left charge conjugation class template. Replacing a fermionic vertex
     * class with this type reverses the first (anti-fermionic Majorana) current, 
     * but doesn't flip the Dirac fermionflow  along the vertex. Its
     * static data is that of the original vertex: */ 
    
    template<class Feynrule_t,bool ferm=Feynrule_t::fermionic>class charge_conj_left
    {
	public:
	    typedef typename Feynrule_t::model_type model_type;
	    static const std::size_t rank=Feynrule_t::rank;
	    static const std::size_t params=Feynrule_t::params;
	    static const bool p_dependent=Feynrule_t::p_dependent;
	    static const bool fermionic=ferm;
	    static const std::size_t tensor_size=Feynrule_t::tensor_size;
	    static const std::size_t sizes[4];
	    static const std::string formula;
    };
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_left<Feynrule_t,ferm>::rank;
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_left<Feynrule_t,ferm>::params;
    template<class Feynrule_t,bool ferm>const bool charge_conj_left<Feynrule_t,ferm>::p_dependent;
    template<class Feynrule_t,bool ferm>const bool charge_conj_left<Feynrule_t,ferm>::fermionic;
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_left<Feynrule_t,ferm>::tensor_size;
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_left<Feynrule_t,ferm>::sizes[4]={Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]};
    template<class Feynrule_t,bool ferm>const std::string charge_conj_left<Feynrule_t,ferm>::formula=Feynrule_t::formula;

    /* For bosonic vertices, the left charge conjugation passes the recursive
     * relations of the original vertex: */

    template<class Feynrule_t>class evaluate< charge_conj_left<Feynrule_t,false> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For fermionic vertices, the left charge conjugation passes the Cc_first,
     * Cc_second, Cc_third and Cc_fourth methods in the original evaluate class,
     * which are the recursive relations of C^* times the Dirac matrix
     * structure: */

    template<class Feynrule_t>class evaluate< charge_conj_left<Feynrule_t,true> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the left charge conjugation evaluation
     * passes the composition of the colour structure and the left charge
     * conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class evaluate<charge_conj_left< compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,charge_conj_left<Feynrule_t2,true> > Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the left charge conjugated vertex
     * cfd-evaluation passes the composition of the colour structure and the
     * left charge conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class cfd_evaluate<charge_conj_left< compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,charge_conj_left<Feynrule_t2,true> > Feynrule_type;
	    typedef typename cfd_evaluate<Feynrule_type>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_type>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_type>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_type>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return cfd_evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::first(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::second(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::third(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta,produced_iters);
	    }
    };

    /* Right charge conjugation class template. Replacing a fermionic vertex
     * class with this type reverses the second (fermionic Majorana) current, 
     * but doesn't flip the Dirac fermion flow along the vertex. Its static data
     * is that of the original vertex: */ 
    

    template<class Feynrule_t,bool ferm=Feynrule_t::fermionic>class charge_conj_right
    {
	public:
	    typedef typename Feynrule_t::model_type model_type;

	    static const std::size_t rank=Feynrule_t::rank;
	    static const std::size_t params=Feynrule_t::params;
	    static const bool p_dependent=Feynrule_t::p_dependent;
	    static const bool fermionic=ferm;
	    static const std::size_t tensor_size=Feynrule_t::tensor_size;
	    static const std::size_t sizes[4];
	    static const std::string formula;
    };
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_right<Feynrule_t,ferm>::rank;
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_right<Feynrule_t,ferm>::params;
    template<class Feynrule_t,bool ferm>const bool charge_conj_right<Feynrule_t,ferm>::p_dependent;
    template<class Feynrule_t,bool ferm>const bool charge_conj_right<Feynrule_t,ferm>::fermionic;
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_right<Feynrule_t,ferm>::tensor_size;
    template<class Feynrule_t,bool ferm>const std::size_t charge_conj_right<Feynrule_t,ferm>::sizes[4]={Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]};
    template<class Feynrule_t,bool ferm>const std::string charge_conj_right<Feynrule_t,ferm>::formula=Feynrule_t::formula;

    /* For bosonic vertices, the right charge conjugation passes the recursive
     * relations of the original vertex: */

    template<class Feynrule_t>class evaluate< charge_conj_right<Feynrule_t,false> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For fermionic vertices, the right charge conjugation passes the first_C,
     * second_C, third_C and fourth_C methods in the original evaluate class,
     * which are the recursive relations of the Dirac matrix structure times C: */

    template<class Feynrule_t>class evaluate< charge_conj_right<Feynrule_t,true> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first_C(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second_C(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third_C(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth_C(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the right charge conjugation evaluation
     * passes the composition of the colour structure and the right charge
     * conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class evaluate<charge_conj_right< compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,charge_conj_right<Feynrule_t2,true> > Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the right charge conjugation
     * cfd-evaluation passes the composition of the colour structure and the
     * right charge conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class cfd_evaluate<charge_conj_right< compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,charge_conj_right<Feynrule_t2,true> > Feynrule_type;
	    typedef typename cfd_evaluate<Feynrule_type>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_type>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_type>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_type>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return cfd_evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::first(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::second(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::third(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta,produced_iters);
	    }
    };

    /* Reversed charge conjugation class template. Replacing a fermionic vertex
     * class with this type flips the Dirac fermion flow along the vertex. Its
     * static data is that of the original vertex: */ 
    
    template<class Feynrule_t,bool ferm=Feynrule_t::fermionic>class reverse
    {
	public:
	    typedef typename Feynrule_t::model_type model_type;

	    static const std::size_t rank=Feynrule_t::rank;
	    static const std::size_t params=Feynrule_t::params;
	    static const bool p_dependent=Feynrule_t::p_dependent;
	    static const bool fermionic=ferm;
	    static const std::size_t tensor_size=Feynrule_t::tensor_size;
	    static const std::size_t sizes[4];
	    static const std::string formula;
    };
    template<class Feynrule_t,bool ferm>const std::size_t reverse<Feynrule_t,ferm>::rank;
    template<class Feynrule_t,bool ferm>const std::size_t reverse<Feynrule_t,ferm>::params;
    template<class Feynrule_t,bool ferm>const bool reverse<Feynrule_t,ferm>::p_dependent;
    template<class Feynrule_t,bool ferm>const bool reverse<Feynrule_t,ferm>::fermionic;
    template<class Feynrule_t,bool ferm>const std::size_t reverse<Feynrule_t,ferm>::tensor_size;
    template<class Feynrule_t,bool ferm>const std::size_t reverse<Feynrule_t,ferm>::sizes[4]={Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]};
    template<class Feynrule_t,bool ferm>const std::string reverse<Feynrule_t,ferm>::formula=Feynrule_t::formula;

    /* For bosonic vertices, the reversal passes the recursive
     * relations of the original vertex: */

    template<class Feynrule_t>class evaluate< reverse<Feynrule_t,false> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For fermionic vertices, the reversal class passes the reversed_first,
     * reversed_second, reversed_third and reversed_fourth methods in the
     * original evaluate class, which are the recursive relations of the
     * reversed Dirac matrix structure: */

    template<class Feynrule_t>class evaluate< reverse<Feynrule_t,true> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the reversal evaluation passes the
     * composition of the colour structure and the reversed inner vertex type:
     * */

    template<class Feynrule_t1,class Feynrule_t2>class evaluate<reverse<compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,reverse<Feynrule_t2,true> > Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the reversal cfd-evaluation passes the
     * composition of the colour structure and the reversed inner vertex type:
     * */

    template<class Feynrule_t1,class Feynrule_t2>class cfd_evaluate<reverse< compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,reverse<Feynrule_t2,true> > Feynrule_type;
	    typedef typename cfd_evaluate<Feynrule_type>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_type>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_type>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_type>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return cfd_evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::first(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::second(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::third(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta,produced_iters);
	    }
    };

    /* Reversed right-charge conjugation class template. Replacing a fermionic
     * vertex class with this type reverses the second (fermionic) current and
     * flips the Dirac fermion flow along the vertex: */ 
    
    template<class Feynrule_t,bool ferm=Feynrule_t::fermionic>class reverse_conj_right
    {
	public:
	    typedef typename Feynrule_t::model_type model_type;

	    static const std::size_t rank=Feynrule_t::rank;
	    static const std::size_t params=Feynrule_t::params;
	    static const bool p_dependent=Feynrule_t::p_dependent;
	    static const bool fermionic=ferm;
	    static const std::size_t tensor_size=Feynrule_t::tensor_size;
	    static const std::size_t sizes[4];
	    static const std::string formula;
    };
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_right<Feynrule_t,ferm>::rank;
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_right<Feynrule_t,ferm>::params;
    template<class Feynrule_t,bool ferm>const bool reverse_conj_right<Feynrule_t,ferm>::p_dependent;
    template<class Feynrule_t,bool ferm>const bool reverse_conj_right<Feynrule_t,ferm>::fermionic;
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_right<Feynrule_t,ferm>::tensor_size;
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_right<Feynrule_t,ferm>::sizes[4]={Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]};
    template<class Feynrule_t,bool ferm>const std::string reverse_conj_right<Feynrule_t,ferm>::formula=Feynrule_t::formula;

    /* For bosonic vertices, the reversed left charge conjugation passes the
     * recursive relations of the original vertex: */

    template<class Feynrule_t>class evaluate< reverse_conj_right<Feynrule_t,false> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For fermionic vertices, the reversed right charge conjugation passes the
     * reversed_first_C, reversed_second_C, reversed_third_C and
     * reversed_fourth_C methods in the original evaluate class, which are the
     * recursive relations of the reversed Dirac matrix structure times C: */

    template<class Feynrule_t>class evaluate< reverse_conj_right<Feynrule_t,true> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_first_C(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_second_C(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_third_C(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::reversed_fourth_C(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the reversed right charge conjugation
     * evaluation passes the composition of the colour structure and the
     * reversed right charge conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class evaluate<reverse_conj_right<compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,reverse_conj_right<Feynrule_t2,true> > Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the reversed right charge conjugation
     * cfd-evaluation passes the composition of the colour structure and the
     * reversed right charge conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class cfd_evaluate<reverse_conj_right< compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,reverse_conj_right<Feynrule_t2,true> > Feynrule_type;
	    typedef typename cfd_evaluate<Feynrule_type>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_type>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_type>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_type>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return cfd_evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::first(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::second(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::third(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta,produced_iters);
	    }
    };

    /* Reversed left-charge conjugation class template. Replacing a fermionic vertex
     * class with this type reverses the first (anti-fermionic) current and flips
     * the Dirac fermion flow along the vertex: */ 
    
    template<class Feynrule_t,bool ferm=Feynrule_t::fermionic>class reverse_conj_left
    {
	public:
	    typedef typename Feynrule_t::model_type model_type;

	    static const std::size_t rank=Feynrule_t::rank;
	    static const std::size_t params=Feynrule_t::params;
	    static const bool p_dependent=Feynrule_t::p_dependent;
	    static const bool fermionic=ferm;
	    static const std::size_t tensor_size=Feynrule_t::tensor_size;
	    static const std::size_t sizes[4];
	    static const std::string formula;
    };
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_left<Feynrule_t,ferm>::rank;
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_left<Feynrule_t,ferm>::params;
    template<class Feynrule_t,bool ferm>const bool reverse_conj_left<Feynrule_t,ferm>::p_dependent;
    template<class Feynrule_t,bool ferm>const bool reverse_conj_left<Feynrule_t,ferm>::fermionic;
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_left<Feynrule_t,ferm>::tensor_size;
    template<class Feynrule_t,bool ferm>const std::size_t reverse_conj_left<Feynrule_t,ferm>::sizes[4]={Feynrule_t::sizes[0],Feynrule_t::sizes[1],Feynrule_t::sizes[2],Feynrule_t::sizes[3]};
    template<class Feynrule_t,bool ferm>const std::string reverse_conj_left<Feynrule_t,ferm>::formula=Feynrule_t::formula;

    /* For bosonic vertices, the reversed right charge conjugation passes the
     * recursive relations of the original vertex: */

    template<class Feynrule_t>class evaluate< reverse_conj_left<Feynrule_t,false> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For fermionic vertices, the left reversed charge conjugation passes the
     * Cc_reversed_first, Cc_reversed_second, Cc_reversed_third and
     * Cc_reversed_fourth methods in the original evaluate class, which are the
     * recursive relations of C* times the reversed Dirac matrix structure: */

    template<class Feynrule_t>class evaluate< reverse_conj_left<Feynrule_t,true> >
    {
	public:
	    typedef Feynrule_t Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_reversed_first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_reversed_second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_reversed_third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::Cc_reversed_fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the reversed left charge conjugation
     * evaluation passes the composition of the colour structure and the
     * reversed left charge conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class evaluate<reverse_conj_left<compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,reverse_conj_left<Feynrule_t2,true> > Feynrule_type;
	    typedef typename evaluate<Feynrule_type>::model_type model_type;
	    typedef typename evaluate<Feynrule_type>::value_type value_type;
	    typedef typename evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename evaluate<Feynrule_type>::size_type size_type;
	    typedef typename evaluate<Feynrule_type>::iterator iterator;
	    typedef typename evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(ARG_LIST)
	    {
		evaluate<Feynrule_type>::first(factor,couplings,iters,momenta);
	    }
	    static void second(ARG_LIST)
	    {
		evaluate<Feynrule_type>::second(factor,couplings,iters,momenta);
	    }
	    static void third(ARG_LIST)
	    {
		evaluate<Feynrule_type>::third(factor,couplings,iters,momenta);
	    }
	    static void fourth(ARG_LIST)
	    {
		evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta);
	    }
    };

    /* For coloured fermionic vertices, the reversed left charge conjugation
     * cfd-evaluation passes the composition of the colour structure and the
     * reversed left charge conjugated inner vertex type: */

    template<class Feynrule_t1,class Feynrule_t2>class cfd_evaluate<reverse_conj_left< compose_vertices<Feynrule_t1,Feynrule_t2>,true> >
    {
	public:
	    typedef compose_vertices<Feynrule_t1,reverse_conj_left<Feynrule_t2,true> > Feynrule_type;
	    typedef typename cfd_evaluate<Feynrule_type>::model_type model_type;
	    typedef typename cfd_evaluate<Feynrule_type>::value_type value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::r_value_type r_value_type;
	    typedef typename cfd_evaluate<Feynrule_type>::tensor_type tensor_type;
	    typedef typename cfd_evaluate<Feynrule_type>::size_type size_type;
	    typedef typename cfd_evaluate<Feynrule_type>::iterator iterator;
	    typedef typename cfd_evaluate<Feynrule_type>::momentum_type momentum_type;

	    static void initialise()
	    {
		cfd_evaluate<Feynrule_type>::initialise();
	    }
	    static std::vector<size_type>& get_index_ranges(size_type n,std::vector<size_type>& v)
	    {
		return cfd_evaluate<Feynrule_type>::get_index_ranges(n,v);
	    }
	    static void first(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::first(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void second(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::second(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void third(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::third(factor,couplings,iters,momenta,produced_iters);
	    }
	    static void fourth(CFD_ARG_LIST)
	    {
		cfd_evaluate<Feynrule_type>::fourth(factor,couplings,iters,momenta,produced_iters);
	    }
    };
}

#include <Camgen/undef_args.h>

#endif /*CAMGEN_CHARGE_CONJ_H_*/

