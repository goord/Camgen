#include "scalar_model.h"
#include <Camgen/CM_algo.h>
#include <Camgen/proc_gen.h>
#include <Camgen/stdrand.h>

using namespace Camgen;

int main(int argc,char *argv[])
{
    bool publisher_mode=false;

    int N_events=1000;
    if(argc>1)
    {
	N_events=std::atoi(argv[1]);
    }

    CM_algorithm<scalar_model,2,4>algo("B,B > A,A,A,A");
    algo.load();
    algo.construct();
    algo.print_tree(std::cout);

    set_channel_init(0,0);
    set_grid_init(0,0);
    set_adaptive_s_sampling(false);
    set_initial_state_type(initial_states::partonic);
    set_phase_space_generator_type(phase_space_generators::uniform);

    set_first_beam_energy(500);
    set_second_beam_energy(500);

    process_generator<scalar_model,2,4,std::random>gen(algo.get_tree_iterator());
    
    gen.generate();
    gen.print(std::cout);
    gen.print_status(std::cout);

    return 0;

    gen.initialise(true);

    gen.print_status(std::cout);

    for(int i=0;i<N_events;++i)
    {
	gen.generate();
	gen.refresh_cross_section();
    }

    gen.print_status(std::cout);

    std::cout<<gen.cross_section()<<std::endl;
}
