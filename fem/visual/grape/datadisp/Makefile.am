# only build this stuff if GRAPE was found
AM_CXXFLAGS = -g -DNDEBUG
if GRAPE
# einfach alle Pakete verwenden, die überhaupt gefunden wurden
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) -lg2c
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

noinst_PROGRAMS = datadisp
datadisp_SOURCES = datadisp.cc 

dist_data_DATA = datadisp.cc printhelp.cc readdata.cc readfile.cc readparams.cc
endif
