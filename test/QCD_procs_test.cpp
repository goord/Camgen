//
// This file is part of the CAMGEN library.
// Copyright (C) 2013 Gijs van den Oord.
// CAMGEN is licensed under the GNU GPL, version 2,
// see COPYING for details.
//

#include <Camgen/QCD.h>
#include <Camgen/CM_algo.h>

/* * * * * * * * * * * * * * * *  * * * * * * * * * *
 * Process-counting test program for QCD amplitudes *
 *                                                  *
 * * * * * * * * ** * * * * * * * * * * * * * * * * */

namespace Camgen
{
    namespace test_utils
    {
	/* QCD process counting facility: */

	class qcd_processes
	{
	    public:

		typedef long long unsigned integer_type;
		typedef std::vector<integer_type>::size_type size_type;

		/* Constructor, initialises the process array up to N final
		 * state partons and f final state flavours: */

		qcd_processes(size_type N,size_type f)
		{
		    initialise(N,f);
		}

		/* Maximal nr. of final state partons supported by the counter
		 * facility: */

		size_type n_max() const
		{
		    return a.size()>0?(a[0].size()):0;
		}

		/* Maximal nr. of final state quark flavours supported by the
		 * process counter: */

		size_type f_max() const
		{
		    return a.size();
		}

		/* Process count for gg initial state: */

		integer_type gg(size_type n,size_type f_out) const
		{
		    return A(n,f_out);
		}

		/* Process count for qqbar initial state: */

		integer_type qqbar(size_type n,size_type f_out,size_type f_in) const
		{
		    return f_in*A(n,f_out);
		}

		/* Process count for gq initial state: */

		integer_type gq(size_type n,size_type f_out,size_type f_in) const
		{
		    return (n==0)?0:(f_in*A(n-1,f_out));
		}

		/* Process count for qq initial state: */

		integer_type qq(size_type n,size_type f_out,size_type f_in) const
		{
		    return (n<2)?0:(f_in*A(n-2,f_out));
		}

		/* Process count for qq' initial state: */

		integer_type qqprime(size_type n,size_type f_out,size_type f_in) const
		{
		    return (n<2)?0:(f_in*(f_in-1)*A(n-2,f_out)/2);
		}

		/* Process count for q'qbar initial state: */

		integer_type qqbarprime(size_type n,size_type f_out,size_type f_in) const
		{
		    return (n<2)?0:(f_in*(f_in-1)*A(n-2,f_out));
		}

		/* Process count for all (ordered) initial parton combinations: */

		integer_type pp(size_type n,size_type f_out,size_type f_in) const
		{
		    return gg(n,f_out)+qqbar(n,f_out,f_in)+2*gq(n,f_out,f_in)+2*qq(n,f_out,f_in)+2*qqprime(n,f_out,f_in)+qqbarprime(n,f_out,f_in);
		}

	    private:

		std::vector< std::vector<integer_type> >a;

		integer_type A(size_type n,size_type f) const
		{
		    return (f<f_max() and n<n_max())?a[f][n]:0;
		}

		void initialise(size_type N,size_type f)
		{
		    if(f>0)
		    {
			initialise(N,f-1);
		    }
		    a.resize(f+1);
		    if(f==0)
		    {
			a[f].resize(N,1);
			return;
		    }
		    else
		    {
			a[f].resize(N);
		    }
		    if(N>0)
		    {
			a[f][0]=1;
		    }
		    if(N>1)
		    {
			a[f][1]=1;
		    }
		    if(N>2)
		    {
			for(size_type n=2;n<N;++n)
			{
			    for(size_type m=0;m<=n/2;++m)
			    {
				a[f][n]+=a[f-1][n-2*m];
			    }
			}
		    }
		}
	};

	/* Function checking whether the number of nonempty pp > nj processes
	 * matches the number calculated by the counter. */

	template<class model_type,std::size_t N>bool check_process_count(const std::string& pin,std::size_t fin,const std::string& pout,std::size_t fout,const qcd_processes* counter=NULL)
	{
	    std::string process = pin + "," + pin + " > ";
	    for(std::size_t n=0;n<N;++n)
	    {
		process.append(pout);
		if(n<N-1)
		{
		    process.append(",");
		}
	    }
	    std::cerr<<"Checking the number of distinct processes for "<<process<<".....";
	    CM_algorithm<model_type,2,N>algo(process);
	    algo.load();
	    algo.construct_trees();
	    algo.remove_empty_processes();
	    qcd_processes::integer_type count;
	    if(counter==NULL)
	    {
		qcd_processes temp_counter(N,fout);
		count=temp_counter.pp(2,fout,fin);
	    }
	    else
	    {
		if(counter->n_max()<N or counter->f_max()<fout)
		{
		    qcd_processes temp_counter(N,fout);
		    count=temp_counter.pp(N,fout,fin);
		}
		else
		{
		    count=counter->pp(N,fout,fin);
		}
	    }
	    qcd_processes::integer_type procs=algo.n_processes();
	    if(count!=procs)
	    {
		std::cerr<<std::endl<<"error: "<<procs<<" processes counted, "<<count<<" processes predicted"<<std::endl;
		return false;
	    }
	    std::cerr<<".....done, "<<count<<" processes counted."<<std::endl;
	    return true;
	}
    }
}

using namespace Camgen;

int main()
{
    license_print::disable();
    std::cout<<"----------------------------------------------------------"<<std::endl;
    std::cout<<"testing process enumeration for pp -> jets in QCD........."<<std::endl;
    std::cout<<"----------------------------------------------------------"<<std::endl;
    
    test_utils::qcd_processes counter(10,10);

    unsigned fj=4;
    unsigned fJ=5;
    unsigned fJt=6;
    unsigned fp=4;
    unsigned fP=5;

    if(!test_utils::check_process_count<QCD,2>("p",fp,"j",fj,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,2>("p",fp,"J",fJ,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,2>("P",fP,"J",fJ,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,2>("p",fp,"Jt",fJt,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,2>("P",fP,"Jt",fJt,&counter))
    {
	return 1;
    }

    if(!test_utils::check_process_count<QCD,3>("p",fp,"j",fj,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,3>("p",fp,"J",fJ,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,3>("P",fP,"J",fJ,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,3>("p",fp,"Jt",fJt,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,3>("P",fP,"Jt",fJt,&counter))
    {
	return 1;
    }

    if(!test_utils::check_process_count<QCD,4>("p",fp,"j",fj,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,4>("p",fp,"J",fJ,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,4>("P",fP,"J",fJ,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,4>("p",fp,"Jt",fJt,&counter))
    {
	return 1;
    }
    if(!test_utils::check_process_count<QCD,4>("P",fP,"Jt",fJt,&counter))
    {
	return 1;
    }
}

