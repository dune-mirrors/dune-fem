# only build this stuff if GRAPE was found
AUTOMAKE_OPTIONS = foreign 1.5

if GRAPE
# einfach alle Pakete verwenden, die überhaupt gefunden wurden

GRIDTYPE=YASPGRID
GRIDDIM= 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) -D$(GRIDTYPE) -DGRIDDIM=$(GRIDDIM)
#noinst_PROGRAMS = gridview
#gridview_SOURCES = gridview.cc 

endif

gridviewdir =  $(includedir)/gridview
gridview_HEADERS = gridview.cc
