#ifndef SCALAR_MODEL_H_
#define SCALAR_MODEL_H_

#include <Camgen/model.h>
#include <Camgen/Minkowski.h>
#include <Camgen/Pauli_basis.h>
#include <Camgen/hel_type.h>

using namespace Camgen;

class scalar_model: public model<scalar_model>
{
    public:

	typedef double value_type;
	typedef Minkowski_type spacetime_type;
	typedef Pauli_basis Dirac_algebra_type;
	typedef helicity_type spin_vector_type;
	
	static const unsigned dimension=4;
	static const int beam_direction=3;

	static const bool coloured=false;
	static const bool continuous_helicities=false;
	static const bool continuous_colours=false;

	static value_type M_h,W_h,M_A,M_B;
	static std::complex<value_type>alpha;

	scalar_model();

	static void set_alpha_s(const double&);
	static void set_QCD_scale(const double&);
};

#endif

