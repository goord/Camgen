//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/ppranknum.h>

#define DELTA(I,i,J,j,n,z)								\
for(size_type b=0;b<N;++b)								\
{											\
    evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);			\
    iters[I]+=sizes[i][I];								\
    iters[J]+=sizes[j][J];								\
}											\
iters[I]-=sizes[i+1][I];								\
iters[J]-=sizes[j+1][J]


#define DDELTA(I,i,J,j,K,k,L,l,n,z)							\
for(size_type a=0;a<N;++a)								\
{											\
    DELTA(K,k,L,l,n,z);									\
    iters[I]+=sizes[i][I];								\
    iters[J]+=sizes[j][J];								\
}											\
iters[I]-=sizes[i+1][I];								\
iters[J]-=sizes[j+1][J]


#define G_to_G_CHAIN(I,J,n,z)								\
for(size_type a=0;a<N*N;++a)								\
{											\
    evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);			\
    iters[I]+=sizes[0][I];								\
    iters[J]+=sizes[0][J];								\
}											\
iters[I]-=sizes[2][I];									\
iters[J]-=sizes[2][J]


#define GG_to_S_CHAIN(I,J,n,z)								\
for(size_type a=0;a<N;++a)								\
{											\
    for(size_type b=0;b<N;++b)								\
    {											\
    	evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);		\
	iters[I]+=sizes[0][I];								\
	iters[J]+=sizes[1][J];								\
    }											\
    iters[J]-=sizes[2][J];								\
    iters[J]+=sizes[0][J];								\
}											\
iters[I]-=sizes[2][I];									\
iters[J]-=sizes[1][J]


#define GG_to_G_CHAIN(I,J,K,n,z)							\
for(size_type a=0;a<N;++a)								\
{											\
    for(size_type b=0;b<N;++b)								\
    {											\
	for(size_type c=0;c<N;++c)							\
	{										\
	    evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);		\
	    iters[I]+=sizes[0][I];							\
	    iters[K]+=sizes[0][K];							\
	}										\
	iters[I]-=sizes[1][I];								\
	iters[J]+=sizes[0][J];								\
    }											\
    iters[I]+=sizes[1][I];								\
    iters[K]-=sizes[2][K];								\
}											\
iters[I]-=sizes[2][I];									\
iters[J]-=sizes[2][J]


#define GGG_to_S_CHAIN(I,J,K,n,z)							\
for(size_type a=0;a<N;++a)								\
{											\
    for(size_type b=0;b<N;++b)								\
    {											\
	for(size_type c=0;c<N;++c)							\
	{										\
	    evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);		\
	    iters[I]+=sizes[1][I];							\
	    iters[K]+=sizes[0][K];							\
	}										\
	iters[I]-=sizes[2][I];								\
	iters[J]+=sizes[0][J];								\
    }											\
    iters[I]+=sizes[0][I];								\
    iters[K]-=sizes[2][K];								\
}											\
iters[I]-=sizes[1][I];									\
iters[J]-=sizes[2][J]


#define GGG_to_G_CHAIN(I,J,K,L,n,z)							\
for(size_type a=0;a<N;++a)								\
{											\
    for(size_type b=0;b<N;++b)								\
    {											\
	for(size_type c=0;c<N;++c)							\
	{										\
	    for(size_type d=0;d<N;++d)							\
	    {										\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);	\
		iters[I]+=sizes[0][I];							\
		iters[L]+=sizes[0][L];							\
	    }										\
	    iters[I]-=sizes[1][I];							\
	    iters[K]+=sizes[0][K];							\
	}										\
	iters[J]+=sizes[0][J];								\
	iters[L]-=sizes[2][L];								\
    }											\
    iters[I]+=sizes[1][I];								\
    iters[K]-=sizes[2][K];								\
}											\
iters[I]-=sizes[2][I];									\
iters[J]-=sizes[2][J]


#define GGGG_to_S_CHAIN(I,J,K,L,n,z)							\
for(size_type a=0;a<N;++a)								\
{											\
    for(size_type b=0;b<N;++b)								\
    {											\
	for(size_type c=0;c<N;++c)							\
	{										\
	    for(size_type d=0;d<N;++d)							\
	    {										\
		evaluate<Feynrule_t>::RANK_NUMERAL(n)(z,couplings,iters,momenta);	\
		iters[I]+=sizes[1][I];							\
		iters[L]+=sizes[0][L];							\
	    }										\
	    iters[I]-=sizes[2][I];							\
	    iters[K]+=sizes[0][K];							\
	}										\
	iters[J]+=sizes[0][J];								\
	iters[L]-=sizes[2][L];								\
    }											\
    iters[I]+=sizes[0][I];								\
    iters[K]-=sizes[2][K];								\
}											\
iters[I]-=sizes[1][I];									\
iters[J]-=sizes[2][J]



