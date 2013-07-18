#include <Camorra/QCD.h>
#include <Camorra/CM_algo.h>
#include <Camorra/uniform_massless_ps.h>
#include <Camorra/QCD_cols.h>
#include <Camorra/rcarry.h>

using namespace Camorra;

int main()
{
    CM_algorithm<QCD,2,10>algo("g,g > g,g,g,g,g,g,g,g,g,g");
    algo.load();
    algo.construct();
    std::cout<<"Nr of Feynman diagrams:"<<algo.count_diagrams()<<std::endl;
    double Ecm=1000;
    algo.set_momentum_generator<uniform_massless_ps,rcarry>(&Ecm);
    algo.set_colour_generator<colour_flow_QCD,rcarry>();
    algo.generate();
    

//    std::ofstream os;
//    os.open("example1.txt");
//    os<<"This is a file containing the gluon colors and momenta that yield the amplitudes quoted for example 1 in the pub's appendix"<<std::endl<<std::endl;
//    os<<"part    cols  \tE                  \tpx                   \tpy                    \tpz                 "<<std::endl;
//    os<<"---------------------------------------------------------------------------------------------------------------------------"<<std::endl;
//    os.precision(10);
//    os.setf(std::ios::scientific,std::ios::floatfield);
//    for(unsigned i=0;i<12;++i)
//    {
//	os<<i<<"\t("<<algo.get_phase_space(i)->colour(0)<<","<<algo.get_phase_space(i)->colour(1)<<")\t"<<algo.get_phase_space(i)->momentum(0);
//	os<<"\t"<<algo.get_phase_space(i)->momentum(1)<<"\t"<<algo.get_phase_space(i)->momentum(2)<<"\t"<<algo.get_phase_space(i)->momentum(3)<<std::endl;
//   }
//    os.close();
    
    
    std::cout<<algo.evaluate()<<std::endl;
    
    algo.get_phase_space(2)->helicity_phase(1)=1.0;			
    algo.get_phase_space(2)->helicity_phase(-1)=1.0;		
    
    std::cout<<algo.evaluate()<<std::endl;			
    
    algo.get_phase_space(3)->helicity_phase(1)=1.0;	
    algo.get_phase_space(3)->helicity_phase(-1)=1.0;
    
    std::cout<<algo.evaluate()<<std::endl;
}
