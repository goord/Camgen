AM_CXXFLAGS =	       -I$(top_srcdir)/include @AM_CXXFLAGS@

lib_LTLIBRARIES =      libCamgen.la libCamROOT.la

libCamgen_la_SOURCES = $(top_srcdir)/src/bipart.cpp		\
		       $(top_srcdir)/src/Dirac_dim.cpp		\
		       $(top_srcdir)/src/EWSM.cpp		\
		       $(top_srcdir)/src/file_utils.cpp		\
		       $(top_srcdir)/src/license_print.cpp	\
		       $(top_srcdir)/src/logstream.cpp		\
		       $(top_srcdir)/src/lower_binom.cpp	\
		       $(top_srcdir)/src/MC_config.cpp		\
		       $(top_srcdir)/src/num_config.cpp		\
		       $(top_srcdir)/src/pdf_wrapper.cpp	\
		       $(top_srcdir)/src/phi3.cpp		\
		       $(top_srcdir)/src/phi34.cpp		\
		       $(top_srcdir)/src/plt_config.cpp		\
		       $(top_srcdir)/src/plt_obj.cpp		\
		       $(top_srcdir)/src/plt_script.cpp		\
		       $(top_srcdir)/src/plt_strm.cpp		\
		       $(top_srcdir)/src/multi_plot.cpp		\
		       $(top_srcdir)/src/QED.cpp		\
		       $(top_srcdir)/src/QCD.cpp		\
		       $(top_srcdir)/src/rcarry.cpp		\
		       $(top_srcdir)/src/scalar_QED.cpp		\
		       $(top_srcdir)/src/SM.cpp			\
		       $(top_srcdir)/src/spin.cpp		\
		       $(top_srcdir)/src/stdrand.cpp		\
		       $(top_srcdir)/src/susy_QED.cpp		\
		       $(top_srcdir)/src/susy_QCD.cpp		\
		       $(top_srcdir)/src/unused.cpp		\
		       $(top_srcdir)/src/WZM.cpp

libCamROOT_la_SOURCES = $(top_srcdir)/src/root_tree.cpp
libCamROOT_la_CXXFLAGS = -I$(top_srcdir)/include $(AM_CXXFLAGS_ROOT)
libCamROOT_la_LDFLAGS = $(AM_LDFLAGS_ROOT)
libCamROOT_la_LIBADD = $(AM_LIBS_ROOT)
