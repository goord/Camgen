#ifndef SMSEESAW_H_
#define SMSEESAW_H_

#include <Camorra/SM_base.h>
#include <Camorra/Minkowski.h>
#include <Camorra/Pauli_basis.h>
#include <Camorra/helicity_type.h>
#include <Camorra/col_flow.h>

using namespace Camorra;

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

