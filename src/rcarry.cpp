//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/rcarry.h>

namespace Camgen
{
    /* Minimal and maximal thrown values: */

    const rcarry::result_type rcarry::min_value;
    const rcarry::result_type rcarry::max_value;

    /* Register of seeds: */

    rcarry::result_type rcarry::seeds[24]={1239423,
                                           222216,
                                           3122644,
                                           523123,
                                           21342,
                                           322845,
                                           82394,
                                           634012,
                                           430458,
                                           31134,
                                           7312345,
                                           8132235,
                                           203496,
                                           8232221,
                                           612342,
                                           123461,
                                           22311267,
                                           371243,
                                           665824,
                                           3211252,
                                           456223,
                                           6123342,
                                           713948,
                                           9223342};

    /* Constructor: */

    rcarry::rcarry():carry(true)
    {
	for(int i=0;i<24;++i)
	{
	    reg[i]=seeds[i];
	}
    }

    /* Throwing operator: */

    rcarry::result_type rcarry::operator()(void)
    {
	result_type y=carry?(reg[9]-reg[23]-1):(reg[9]-reg[23]);
	carry=(y<0);
	result_type x=carry?(y+max_value):y;
	for(int i=23;i!=0;--i)
	{
	    reg[i]=reg[i-1];
	}
	reg[0]=x;
	return x;
    }
}

