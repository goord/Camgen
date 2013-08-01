#ifndef SMSEESAW_H_
#define SMSEESAW_H_

#include <Camgen/SM_base.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/hel_type.h>
#include <Camgen/col_flow.h>

using namespace Camgen;

class SMseesaw: public SM_base<SMseesaw,double>
{
    public:
	typedef double value_type;
	typedef Minkowski_type spacetime_type;
	typedef Pauli_basis Dirac_algebra_type;
	typedef helicity_type spin_vector_type;
	typedef colour_flow colour_treatment;
	
	static const unsigned dimension=4;
	static const unsigned N_c=3;
	static const int beam_direction=3;

	static const bool coloured=true;
	static const bool continuous_helicities=false;
	static const bool continuous_colours=false;

	static value_type M_N;

	SMseesaw();
};

#endif

