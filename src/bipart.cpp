//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <algorithm>
#include <Camgen/bipart.h>

namespace Camgen
{
    /* Constructor: */

    bi_partition::bi_partition(const std::set<index_type>& S):n(S.size()),kappa(n,0),maxima(n,0)
    {
	size_type part_size=(n+1)/2;
	if(n!=0)
	{
	    for(size_type i=n-part_size+1;i<n;++i)
	    {
		kappa[i]=i+part_size-n;
		maxima[i]=i+part_size-n;
	    }
	    do
	    {
		add_partition(S);
	    }
	    while(next_index_set(part_size));
	}
    }

    /* Partition construction iteration: */
    
    void bi_partition::add_partition(const std::set<size_type>& S)
    {
	size_type i=0;
	std::vector< std::pair<std::size_t,std::size_t> >part(maxima[n-1]-maxima[0]+1);
	std::vector<int>fillings(maxima[n-1]-maxima[0]+1,0);
	for(std::set<std::size_t>::const_iterator it=S.begin();it!=S.end();++it)
	{
	    size_type m=std::count(kappa.begin(),kappa.end(),kappa[i]);
	    if(m<1 or m>2)
	    {
		return;
	    }
	    if(fillings[kappa[i]]==0)
	    {
		part[kappa[i]].first=*it;
	    }
	    else
	    {
		part[kappa[i]].second=*it;
	    }
	    ++fillings[kappa[i]];
	    ++i;
	}
	for(size_type i=0;i<fillings.size();++i)
	{
	    if(fillings[i]!=2)
	    {
		part[i].second=part[i].first;
	    }
	}
	this->push_back(part);
    }

    /* Partition construction iteration helper: */

    bool bi_partition::next_index_set(size_type p)
    {
	for(size_type i=n-1;i!=0;--i)
	{
	    if(kappa[i]<=maxima[i-1] and kappa[i]<p-1)
	    {
		++(kappa[i]);
		maxima[i]=std::max(maxima[i],kappa[i]);
		for(size_type j=i+1;j!=n-p+maxima[i]+1;++j)
		{
		    kappa[j]=kappa[0];
		    maxima[j]=maxima[i];
		}
		for(size_type j=n-p+maxima[i]+1;j!=n;++j)
		{
		    kappa[j]=p-n+j;
		    maxima[j]=p-n+j;
		}
		return true;
	    }
	}
	return false;
    }

    /* Partition printing facility: */

    std::ostream& operator << (std::ostream& os,const bi_partition& P)
    {
	for(bi_partition::size_type i=0;i<P.size();++i)
	{
	    os<<'[';
	    for(bi_partition::size_type j=0;j<P[i].size()-1;++j)
	    {
		os<<'('<<P[i][j].first<<','<<P[i][j].second<<"),";
	    }
	    os<<'('<<P[i].back().first<<','<<P[i].back().second<<")]"<<std::endl;
	}
	return os;
    }
}

