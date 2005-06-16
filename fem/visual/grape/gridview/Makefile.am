# only build this stuff if GRAPE was found
AUTOMAKE_OPTIONS = foreign 1.5

if GRAPE
# einfach alle Pakete verwenden, die überhaupt gefunden wurden
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

noinst_PROGRAMS = gridview
gridview_SOURCES = gridview.cc
endif
