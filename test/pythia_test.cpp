//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <config.h>
#include <Camgen/file_utils.h>
#include <Camgen/SM.h>
#include <Camgen/stdrand.h>
#include <Camgen/evtgen_fac.h>

#if HAVE_PYTHIA8_H
#include <Camgen/Pythia_if.h>
#include <Pythia8/Pythia.h>
#endif

int main()
{
    Camgen::license_print::disable();

#if HAVE_PYTHIA8_H

    typedef Camgen::SM model_type;
    typedef std::random rn_engine;
    typedef std::size_t size_type;
    typedef double value_type;

    size_type nevts=100;
    Camgen::file_utils::create_directory("test_output/pythia_test");

    Camgen::log.enable_level=Camgen::log_level::error;
    Camgen::set_initial_state_type(Camgen::initial_states::proton_proton);
    Camgen::set_phase_space_generator_type(Camgen::phase_space_generators::recursive);
    std::string process("p,p > l-,nubar,h0");
    std::string fname("test_output/pythia_test/pp_lnh.lhe");
    value_type E1=5000;
    value_type E2=5000;
    std::cerr<<"Checking pythia interface for "<<process<<"............"<<std::endl;
    std::cerr.flush();
    Camgen::CM_algorithm<model_type,2,3>algo(process);
    algo.load();
    algo.construct_trees();
    Camgen::set_beam_energy(-1,E1);
    Camgen::set_beam_energy(-2,E2);
    Camgen::event_generator_factory<model_type,2,3,rn_engine> factory;
    Camgen::event_generator<model_type,2,3,rn_engine>* evt_gen=factory.create_generator(algo);
    evt_gen->refresh_Ecm();
    pre_initialise(evt_gen,true);

    Camgen::Pythia_interface<model_type,3> pif(evt_gen,1,true,1);
    Pythia8::Pythia pythia;
    pythia.readString("Beams:frameType = 5");
    pythia.setLHAupPtr(static_cast<Pythia8::LHAup*>(&pif));
    Pythia8::LHAupFromPYTHIA8 pyLHA(&pythia.process,&pythia.info);
    pyLHA.openLHEF(fname);
    pythia.init();
    pyLHA.setInit();
    pyLHA.initLHEF();
    for(size_type i=0;i<nevts;++i)
    {
        if(!pythia.next())
        {
            return 1;
        }
        pyLHA.setEvent();
        pyLHA.eventLHEF();
    }
    pyLHA.closeLHEF(true);
    delete evt_gen;
    std::cerr<<"done"<<std::endl;
    return 0;
    Camgen::log.enable_level=Camgen::log_level::warning;
    
#else

    std::cerr<<"Skipping Pyhia-8 checks: no Pythia-8 install location provided."<<std::endl;
    return 0;

#endif /* HAVE_PYTHIA8_H*/
}
