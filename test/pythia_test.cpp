//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <config.h>
#include <Camgen/Pythia_if.h>
#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>

#if HAVE_PYTHIA8_H
#include <Pythia8/Pythia.h>
#endif

int main()
{
    Camgen::license_print::disable();

#if HAVE_PYTHIA8_H

    typedef Camgen::SM model_type;
    typedef std::random rn_engine;
//    typedef std::size_t size_type;
    typedef double value_type;

    Camgen::log.enable_level=Camgen::log_level::error;
    Camgen::set_initial_state_type(Camgen::initial_states::proton_proton);
    Camgen::set_phase_space_generator_type(Camgen::phase_space_generators::recursive);
    std::string process("p,p > e-,nu_ebar,j");
    std::string fname("test_output/ascii_test/pp_2l2n");
    value_type E1=1000;
    value_type E2=1000;
    std::cerr<<"Checking pythia interface for "<<process<<"............";
    std::cerr.flush();
    Camgen::CM_algorithm<model_type,2,3>algo(process);
    algo.load();
    algo.construct_trees();
    Camgen::set_beam_energy(-1,E1);
    Camgen::set_beam_energy(-2,E2);
    Camgen::event_generator_factory<model_type,2,3,rn_engine> factory;
    Camgen::event_generator<model_type,2,3,rn_engine>* evt_gen=factory.create_generator(algo);

    Camgen::Pythia_interface<model_type,3> pif(evt_gen,1,true,1);
    Pythia8::Pythia parton_shower;
    parton_shower.setLHAupPtr(static_cast<Pythia8::LHAup*>(&pif));
    parton_shower.init();
    parton_shower.next();

    delete evt_gen;
    Camgen::log.enable_level=Camgen::log_level::warning;

    
#else

    std::cerr<<"Skipping Pyhia8 checks: no Pythia8 install location proided."<<std::endl;

#endif /* HAVE_PYTHIA8_H*/
}
