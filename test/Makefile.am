AM_CPPFLAGS =		    	-O3 -pipe -I$(top_srcdir)/include

EXTRA_DIST =		    	plots

if USE_LHAPDF
AM_LDFLAGS =			-L$(LHAPDF_DIR)/lib -lLHAPDF
endif

check_LTLIBRARIES =	    	libCamtest.la

CLEANFILES =                	*.o *.dat *.plt *.log

check_PROGRAMS =		bitstring_test		\
				Lorentz_test		\
		 	    	phi3_graphs_test	\
		 	    	phi34_graphs_test	\
				QCD_procs_test		\
		 	    	Pauli_basis_test	\
		 	    	Weyl_basis_test		\
		 	    	Pauli_spinor_test	\
		 	    	Weyl_spinor_test	\
		 	    	m_spinor_test		\
		 	    	discr_cols_test		\
		 	    	cont_cols_test		\
		 	    	QEDch_test		\
		 	    	QEDdh_test		\
		 	    	QCDcc_test		\
		 	    	QCDdc_test		\
		 	    	Parke_Taylor_cc_test	\
		 	    	Parke_Taylor_dc_test	\
		 	    	SM_test			\
		 		SM_params_test		\
		 		susy_QED_test		\
		 		susy_QCDcc_test		\
		 		susy_QCDdc_test		\
		 		susy_Kunszt_test	\
		 		MC_hel_test		\
		 		MC_col_test		\
		 		MC_gen_test		\
		 		s_int_test		\
				rambo_test		\
		 		ps_tree_test		\
		 		ps_tree_reverse_test	\
		 		ps_tree_shat_test	\
		 		ps_tree_reverse_shat_test\
				ps_tree_adapt_test	\
				ps_tree_grid_test	\
		 		parni_test		\
		 		LHAPDF_test		\
				ps_vars_test		\
				save_load_test		\
		 		speed_test		

check_HEADERS =	    	    	ps_gen_tester.h		\
				QEDPbdh.h		\
			    	QEDPbch.h		\
			    	QEDWbdh.h		\
			    	QEDWbch.h		\
			    	QEDPbKsdh.h		\
			    	QEDPbpsdh.h		\
			    	QEDWbKsdh.h		\
			    	QEDWbpsdh.h		\
			    	QCDPbchabcc.h		\
			    	QCDPbchabcc_clone.h	\
			    	QCDPbchabdc.h		\
			    	QCDPbchcfcc.h		\
			    	QCDPbchcfcc_clone.h	\
			    	QCDPbchcfdc.h		\
			    	QCDPbchcfdc_clone.h	\
			    	QCDPbdhcfcc.h		\
			    	QCDPbdhcfdc.h		\
			    	SMWbhsch.h		\
			    	SMPbKsch.h		\
				susy_QEDWb.h		\
				susy_QCDPbchabcc.h	\
				susy_QCDPbchcfcc.h	\
				susy_QCDPbchcfdc.h	\
				susy_QCDPbdhcfcc.h	\
				s_int_tester.h		\
				test_gen.h

libCamtest_la_SOURCES =		QEDPbdh.cpp		\
				QEDPbch.cpp		\
			    	QEDWbdh.cpp		\
			    	QEDWbch.cpp		\
			    	QEDPbKsdh.cpp		\
			    	QEDPbpsdh.cpp		\
			    	QEDWbKsdh.cpp		\
			    	QEDWbpsdh.cpp		\
				QCDPbchabcc.cpp		\
			    	QCDPbchabcc_clone.cpp	\
			    	QCDPbchabdc.cpp		\
			    	QCDPbchcfcc.cpp		\
			    	QCDPbchcfcc_clone.cpp	\
			    	QCDPbchcfdc.cpp		\
			    	QCDPbchcfdc_clone.cpp	\
			    	QCDPbdhcfcc.cpp		\
			    	QCDPbdhcfdc.cpp		\
			    	SMWbhsch.cpp		\
			    	SMPbKsch.cpp		\
				susy_QEDWb.cpp		\
				susy_QCDPbchabcc.cpp	\
				susy_QCDPbchcfcc.cpp	\
				susy_QCDPbchcfdc.cpp	\
				susy_QCDPbdhcfcc.cpp

bitstring_test_SOURCES =	bitstring_test.cpp
bitstring_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

Lorentz_test_SOURCES =		Lorentz_test.cpp
Lorentz_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

phi3_graphs_test_SOURCES = 	phi3_graphs_test.cpp
phi3_graphs_test_LDADD =	$(top_srcdir)/lib/libCamgen.la

phi34_graphs_test_SOURCES = 	phi34_graphs_test.cpp
phi34_graphs_test_LDADD =	$(top_srcdir)/lib/libCamgen.la

QCD_procs_test_SOURCES =	QCD_procs_test.cpp
QCD_procs_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

Pauli_basis_test_SOURCES = 	Pauli_basis_test.cpp
Pauli_basis_test_LDADD =	$(top_srcdir)/lib/libCamgen.la

Weyl_basis_test_SOURCES =	Weyl_basis_test.cpp
Weyl_basis_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

Pauli_spinor_test_SOURCES =	Pauli_spinor_test.cpp
Pauli_spinor_test_LDADD =	$(top_srcdir)/lib/libCamgen.la

Weyl_spinor_test_SOURCES =	Weyl_spinor_test.cpp
Weyl_spinor_test_LDADD =	$(top_srcdir)/lib/libCamgen.la

m_spinor_test_SOURCES =		m_spinor_test.cpp
m_spinor_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

discr_cols_test_SOURCES = 	discr_cols_test.cpp
discr_cols_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

cont_cols_test_SOURCES = 	cont_cols_test.cpp
cont_cols_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

QEDch_test_SOURCES = 		QEDch_test.cpp
QEDch_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

QEDdh_test_SOURCES = 		QEDdh_test.cpp
QEDdh_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

QCDcc_test_SOURCES = 		QCDcc_test.cpp
QCDcc_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

QCDdc_test_SOURCES = 		QCDdc_test.cpp
QCDdc_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

Parke_Taylor_cc_test_SOURCES = 	Parke_Taylor_cc_test.cpp
Parke_Taylor_cc_test_LDADD =	$(top_srcdir)/lib/libCamgen.la libCamtest.la

Parke_Taylor_dc_test_SOURCES = 	Parke_Taylor_dc_test.cpp
Parke_Taylor_dc_test_LDADD =	$(top_srcdir)/lib/libCamgen.la libCamtest.la

SM_test_SOURCES	=		SM_test.cpp
SM_test_LDADD =			$(top_srcdir)/lib/libCamgen.la libCamtest.la

SM_params_test_SOURCES =	SM_params_test.cpp
SM_params_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

susy_QED_test_SOURCES =		susy_QED_test.cpp
susy_QED_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

susy_QCDcc_test_SOURCES =	susy_QCDcc_test.cpp
susy_QCDcc_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

susy_QCDdc_test_SOURCES =	susy_QCDdc_test.cpp
susy_QCDdc_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

susy_Kunszt_test_SOURCES =	susy_Kunszt_test.cpp
susy_Kunszt_test_LDADD =	$(top_srcdir)/lib/libCamgen.la libCamtest.la

MC_hel_test_SOURCES = 		MC_hel_test.cpp
MC_hel_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

MC_col_test_SOURCES = 		MC_col_test.cpp
MC_col_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

MC_gen_test_SOURCES =   	MC_gen_test.cpp
MC_gen_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

s_int_test_SOURCES =   	    	s_int_test.cpp
s_int_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

rambo_test_SOURCES =	    	rambo_test.cpp
rambo_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

ps_tree_test_SOURCES =	    	ps_tree_test.cpp
ps_tree_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la

ps_tree_reverse_test_SOURCES =	ps_tree_reverse_test.cpp
ps_tree_reverse_test_LDADD =	$(top_srcdir)/lib/libCamgen.la libCamtest.la

ps_tree_shat_test_SOURCES =	ps_tree_shat_test.cpp
ps_tree_shat_test_LDADD =	$(top_srcdir)/lib/libCamgen.la

ps_tree_reverse_shat_test_SOURCES = ps_tree_reverse_shat_test.cpp
ps_tree_reverse_shat_test_LDADD = $(top_srcdir)/lib/libCamgen.la

ps_tree_adapt_test_SOURCES = 	ps_tree_adapt_test.cpp
ps_tree_adapt_test_LDADD = 	$(top_srcdir)/lib/libCamgen.la

ps_tree_grid_test_SOURCES = 	ps_tree_grid_test.cpp
ps_tree_grid_test_LDADD = 	$(top_srcdir)/lib/libCamgen.la

parni_test_SOURCES =		parni_test.cpp
parni_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

LHAPDF_test_SOURCES =		LHAPDF_test.cpp
LHAPDF_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

ps_vars_test_SOURCES =		ps_vars_test.cpp
ps_vars_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

save_load_test_SOURCES =	save_load_test.cpp
save_load_test_LDADD =		$(top_srcdir)/lib/libCamgen.la

speed_test_SOURCES =		speed_test.cpp
speed_test_LDADD =		$(top_srcdir)/lib/libCamgen.la libCamtest.la


TESTS =			    	speed_test		\
			    	Parke_Taylor_cc_test	\
				bitstring_test		\
				Lorentz_test		\
			    	phi3_graphs_test	\
			    	phi34_graphs_test	\
				QCD_procs_test		\
			    	Pauli_basis_test	\
			    	Weyl_basis_test		\
			    	m_spinor_test		\
			    	Pauli_spinor_test	\
			    	Weyl_spinor_test	\
			    	discr_cols_test		\
			    	cont_cols_test		\
			    	QEDch_test		\
			    	QEDdh_test		\
			    	QCDcc_test		\
			    	QCDdc_test		\
			    	Parke_Taylor_dc_test	\
			    	SM_test			\
				SM_params_test		\
				susy_QED_test		\
				susy_QCDcc_test		\
				susy_QCDdc_test		\
				susy_Kunszt_test	\
				MC_hel_test		\
				MC_col_test		\
				MC_gen_test		\
			    	s_int_test		\
			    	rambo_test		\
			    	ps_tree_test		\
				ps_tree_reverse_test	\
				ps_tree_shat_test	\
				ps_tree_reverse_shat_test\
				ps_tree_adapt_test	\
				ps_tree_grid_test	\
				parni_test		\
				LHAPDF_test		\
				ps_vars_test		\
				save_load_test

