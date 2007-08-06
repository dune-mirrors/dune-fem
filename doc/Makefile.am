# $Id: Makefile.am 4908 2007-04-18 09:12:01Z christi $

# also build these sub directories
SUBDIRS = doxygen

# automatically create these web pages
all: $(PAGES)

# setting like in dune-web
CURDIR=doc
# position of the web base directory, 
# relative to $(CURDIR)
BASEDIR=..

# install the html pages
docdir=$(datadir)/doc/dune-fem

# include rules for wml -> html transformation
include $(top_srcdir)/am/webstuff

# remove html pages on ``make  clean''
CLEANFILES = $(PAGES)

# include further rules needed by Dune
include $(top_srcdir)/am/global-rules
