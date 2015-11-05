SUBDIRS = dune lib m4 doc

# we need the module file to be able to build via dunecontrol
EXTRA_DIST = dune.module modules.txt

DISTCHECK_CONFIGURE_FLAGS = --with-dune="$(DUNEROOT)" CXX="$(CXX)" CC="$(CC)"

SOURCEFILE_MASK = -name \*.h -o -name \*.hh -o -name \*.c -o -name \*.cc

# add library name, authors list, and license to progam files
dist-hook:
	if test -f scripts/addauthors.sh ; then \
		find $(distdir) \( $(SOURCEFILE_MASK) \) -exec scripts/addauthors.sh ./AUTHORS {} \; ; \
	fi

include $(top_srcdir)/am/top-rules
include $(top_srcdir)/am/global-rules
