#include "SMseesaw.h"

using namespace Camgen;

const unsigned SMseesaw::dimension;
const unsigned SMseesaw::N_c;
const int SMseesaw::beam_direction;
const bool SMseesaw::coloured;
const bool SMseesaw::continuous_helicities;
const bool SMseesaw::continuous_colours;

SMseesaw::value_type SMseesaw::M_N=1000;

SMseesaw::SMseesaw()
{
    add_fermion("~nu_e",&M_N);

    add_vertex<vffR>("W+","~nu_e","e-",&SM_base<SMseesaw,double>::Wnee);
    add_vertex<vffR>("W-","e+","~nu_e",&SM_base<SMseesaw,double>::Wnee);
}


