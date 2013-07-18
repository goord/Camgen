#include "SMseesaw.h"
#include <Camorra/CM_algo.h>
#include <Camorra/uniform_massless_ps.h>
#include <Camorra/uniform_hels.h>
#include <Camorra/QCD_cols.h>
#include <Camorra/rcarry.h>

using namespace Camorra;

int main()
{
    CM_algorithm<SMseesaw,2,4>algo("d,ubar > u,dbar,e-,e-",5);
    algo.load();
    algo.construct();
    algo.print_tree(std::cout);
    std::cout<<"Nr of Feynman diagrams:"<<algo.count_diagrams()<<std::endl;
    double Ecm=500;
    algo.set_momentum_generator<uniform_massless_ps,rcarry>(&Ecm);
    algo.set_colour_generator<colour_flow_QCD,rcarry>();
    algo.generate();

//    std::ofstream os;
//    os.open("example2.txt");
//    os<<"This is a file containing the colors and momenta that yield the amplitude quoted for example 2 in the pub's appendix"<<std::endl<<std::endl;
//    os<<"part    col\tE                  \tpx                   \tpy                    \tpz                 "<<std::endl;
//    os<<"--------------------------------------------------------------------------------------------------------------------"<<std::endl;
//    os.precision(10);
//    os.setf(std::ios::scientific,std::ios::floatfield);
//    for(unsigned i=0;i<4;++i)
//    {
//	os<<i<<"\t("<<algo.get_phase_space(i)->colour(0)<<")\t"<<algo.get_phase_space(i)->momentum(0);
//	os<<"\t"<<algo.get_phase_space(i)->momentum(1)<<"\t"<<algo.get_phase_space(i)->momentum(2)<<"\t"<<algo.get_phase_space(i)->momentum(3)<<std::endl;
//    }
//    for(unsigned i=4;i<6;++i)
//    {
//	os<<i<<"\t( )\t"<<algo.get_phase_space(i)->momentum(0);
//	os<<"\t"<<algo.get_phase_space(i)->momentum(1)<<"\t"<<algo.get_phase_space(i)->momentum(2)<<"\t"<<algo.get_phase_space(i)->momentum(3)<<std::endl;
//    }
//    os.close();
    
    algo.sum_spins();
    std::cout<<algo.evaluate_spin_sum()<<std::endl;
}
