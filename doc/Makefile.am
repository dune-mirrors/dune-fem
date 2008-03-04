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

if ! BUILD_DOCS
# add tag to notify that dist has been build without documentation 
dist-hook:
	echo "# No documentation included in distribution! " > $(distdir)/$(DOCUMENTATION_TAG_FILE) 
endif


# include further rules needed by Dune
include $(top_srcdir)/am/global-rules
