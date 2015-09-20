#include "scalar_model.h"
#include <Camgen/sss.h>
#include <Camgen/ssss.h>

using namespace Camgen;

const unsigned scalar_model::dimension;
const int scalar_model::beam_direction;
const bool scalar_model::coloured;
const bool scalar_model::continuous_helicities;
const bool scalar_model::continuous_colours;

scalar_model::value_type scalar_model::M_h=125;
scalar_model::value_type scalar_model::W_h=5;
scalar_model::value_type scalar_model::M_A=5;
scalar_model::value_type scalar_model::M_B=10;
std::complex<scalar_model::value_type> scalar_model::alpha(0,0.5);

scalar_model::scalar_model()
{
    add_scalar("h",&M_h,&W_h);
    add_scalar("A",&M_A);
    add_scalar("B",&M_B);

    add_vertex<sss>("h","A","A",&alpha);
    add_vertex<ssss>("h","h","B","B",&alpha);
}

void scalar_model::set_alpha_s(const double& s)
{
    return;
}

void scalar_model::set_QCD_scale(const double& s)
{
    return;
}
