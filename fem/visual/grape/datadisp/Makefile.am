# only build this stuff if GRAPE was found
#CXXFLAGS = -g -DNDEBUG
if GRAPE
# einfach alle Pakete verwenden, die überhaupt gefunden wurden
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

noinst_PROGRAMS = datadisp paramdisp
datadisp_SOURCES = datadisp.cc 
paramdisp_SOURCES = paramdisp.cc 
endif
