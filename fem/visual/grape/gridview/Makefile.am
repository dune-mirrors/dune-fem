# only build this stuff if GRAPE was found
AUTOMAKE_OPTIONS = foreign 1.5

if GRAPE
# einfach alle Pakete verwenden, die überhaupt gefunden wurden
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

noinst_PROGRAMS = gridview gridtest
#noinst_PROGRAMS = gridtest
gridview_SOURCES = gridview.cc 
gridtest_SOURCES = gridtest.cc 
endif
