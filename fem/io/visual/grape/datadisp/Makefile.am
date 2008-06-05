# only build this stuff if GRAPE was found
# einfach alle Pakete verwenden, die überhaupt gefunden wurden
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

datadispdir =  $(includedir)/dune/fem/io/visual/grape/datadisp
datadisp_HEADERS = dataconvert.cc \
		   datadisp.cc \
		   errordisplay.hh \
		   grcommon.hh \
                   printhelp.cc \
		   programtemplate.cc \
		   readiotupledata.cc \
                   readioparams.cc
