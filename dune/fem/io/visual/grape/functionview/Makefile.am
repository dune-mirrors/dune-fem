# only build this stuff if GRAPE was found
AUTOMAKE_OPTIONS = foreign 1.5

# einfach alle Pakete verwenden, die überhaupt gefunden wurden
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)

#noinst_PROGRAMS =  functionview
#functionview_SOURCES = functionview.cc 

#EXTRA_DIST = functionview.cc

include $(top_srcdir)/am/global-rules

