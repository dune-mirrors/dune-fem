# $Id: Makefile.am 4561 2006-05-17 14:47:35Z mblatt $

if BUILD_DOCS
PAGES = codingstyle.html portability.html \
  tests.html features.html indexidsets.html
endif

all: $(PAGES)

# place in website
BASEDIR=../..
CURDIR=doc/devel

# install the html pages
docdir=$(datadir)/doc/dune-common/devel
doc_DATA = $(PAGES)
EXTRA_DIST= $(PAGES)

include $(top_srcdir)/am/webstuff

CLEANFILES = $(PAGES)

include $(top_srcdir)/am/global-rules
