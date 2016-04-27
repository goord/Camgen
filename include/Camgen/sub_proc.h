//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

/*! \file sub_proc.h
    \brief Class containing sub-process particles.
 */

#ifndef SUB_PROC_H_
#define SUB_PROC_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Class holding incoming and outgoing sub-process particles.  *
 *                                                             *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <Camgen/particle.h>
#include <Camgen/vector.h>

namespace Camgen
{
    /// Class holding in- and outgoing particle pointers.
    
    template<class model_t,std::size_t N_in,std::size_t N_out>class sub_process
    {
        public:

            /* Type definitions: */

            typedef model_t model_type;
            typedef typename model_type::value_type value_type;
            typedef particle<model_type> particle_type;
            typedef typename vector<const particle_type*,N_out>::size_type size_type;

            /// Constructor.
            
            sub_process(const vector<const particle_type*,N_in>& p_in,const vector<const particle_type*,N_out>& p_out)
            {
                for(size_type i=0;i<N_in;++i)
                {
                    incoming_particles[i]=p_in[i];
                }
                for(size_type i=0;i<N_out;++i)
                {
                    outgoing_particles[i]=p_out[i];
                }
            }

            /// Cloning method.

            virtual sub_process<model_t,N_in,N_out>* clone() const
            {
                return new sub_process<model_t,N_in,N_out>(incoming_particles,outgoing_particles);
            }

            /// Returns the i-th incoming particle.

            const particle_type* get_particle_in(size_type i) const
            {
                return i<N_in?incoming_particles[i]:NULL;
            }

            /// Returns the i-th outgoing particle.

            const particle_type* get_particle_out(size_type i) const
            {
                return i<N_out?outgoing_particles[i]:NULL;
            }

            /// Returns the i-th particle.

            const particle_type* get_particle(int i) const
            {
		return (i<0)?get_particle_in(-i-1):get_particle_out(i-1);
            }

            /// Returns the i-th incoming particle pdg id.

            int id_in(size_type i) const
            {
                const particle<model_t>* p=get_particle_in(i);
                return p==NULL?0:p->pdg_id();
            }

            /// Returns the i-th outgoing particle pdg id.

            int id_out(size_type i) const
            {
                const particle<model_t>* p=get_particle_out(i);
                return p==NULL?0:p->pdg_id();
            }

            /// Returns the i-th particle pdg id.

            int id(int i) const 
            {
		return (i<0)?id_in(-i-1):id_out(i-1);
            }

            /// Returns the i-th incoming particle name.

            const std::string& name_in(size_type i) const
            {
                const particle<model_t>* p=get_particle_in(i);
                return p==NULL?"":p->get_name();
            }

            /// Returns the i-th outgoing particle name.

            const std::string& name_out(size_type i) const
            {
                const particle<model_t>* p=get_particle_out(i);
                return p==NULL?"":p->get_name();
            }

            /// Returns the i-th particle name.

            const std::string& name(int i) const 
            {
		return (i<0)?name_in(-i-1):name_out(i-1);
            }

            /// Returns the i-th incoming particle mass.

            value_type m_in(size_type i) const
            {
                const particle<model_t>* p=get_particle_in(i);
                return p==NULL?"":p->get_mass();
            }

            /// Returns the i-th outgoing particle mass.

            value_type m_out(size_type i) const
            {
                const particle<model_t>* p=get_particle_out(i);
                return p==NULL?"":p->get_mass();
            }

            /// Returns the i-th particle mass.

            value_type m(int i) const 
            {
		return (i<0)?m_in(-i-1):m_out(i-1);
            }

        private:

            vector<const particle_type*,N_in> incoming_particles;
            vector<const particle_type*,N_out> outgoing_particles;
    };
}

#endif /*SUB_PROC_H_*/
