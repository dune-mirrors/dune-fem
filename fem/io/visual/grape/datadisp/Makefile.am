# only build this stuff if GRAPE was found
if GRAPE
# einfach alle Pakete verwenden, die überhaupt gefunden wurden
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
endif

datadispdir =  $(includedir)/datadisp
datadisp_HEADERS = leng.hh       readdata.cc  readparams.cc\
readtupparams.cc datadisp.cc  printhelp.cc  readfile.cc  readtupledata.cc

#noinst_PROGRAMS = datadisp
#datadisp_SOURCES = datadisp.cc 
