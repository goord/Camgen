AM_CPPFLAGS =		    	-O3 -pipe -I$(top_srcdir)/include
if USE_LHAPDF
AM_LDFLAGS = 			-L$(top_srcdir)/lib -lCamgen -L$(LHAPDF_DIR)/lib -lLHAPDF
else
AM_LDFLAGS = 			-L$(top_srcdir)/lib -lCamgen
endif

noinst_PROGRAMS = 		example1	\
				example2

example1_SOURCES =		example1.cpp

example2_SOURCES =		example2.cpp	\
				SMseesaw.cpp
