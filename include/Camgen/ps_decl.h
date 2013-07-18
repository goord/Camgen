//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#ifndef CAMGEN_PS_DECL_H_
#define CAMGEN_PS_DECL_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Forward declarations and static data for phase space generation algorithms. *
 *                                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace Camgen
{
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class momentum_channel;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class momentum_channel_comp;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class particle_channel;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class particle_channel_comp;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_branching;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t>class s_branching;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t>class backward_s_branching;
    template<class model_t,std::size_t N_out,class rng_t,class spacetime_t>class t_branching;
    template<class model_t,std::size_t N_out,class rng_t,class spacetime_t>class backward_t_branching;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class sum_branching;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class backward_sum_branching;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t>class ps_tree;
    template<class model_t,std::size_t N_in,std::size_t N_out,class rng_t,class spacetime_t>class ps_factory;
}

#endif /*CAMGEN_PS_DECL_H_*/

