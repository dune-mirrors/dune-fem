# $Id$

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# don't follow the full GNU-standard; we definitely need automake >= 1.5
#AUTOMAKE_OPTIONS = foreign 1.5

# do everything in the subdirs
SUBDIRS = doc fem lib

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dune-fem.pc

DISTCHECK_CONFIGURE_FLAGS = --with-dune-common="$(DUNE_COMMON_ROOT)" --with-dune-grid="$(DUNE_GRID_ROOT)" CXX="$(CXX)" CC="$(CC)" DOXYGEN="true"

include $(top_srcdir)/am/global-rules

# add library name, authors list, and license to progam files
dist-hook:
	find $(distdir) -name "*.h" -exec ./addauthors.sh ./AUTHORS {} \;
	find $(distdir) -name "*.hh" -exec ./addauthors.sh ./AUTHORS {} \;
	find $(distdir) -name "*.c" -exec ./addauthors.sh ./AUTHORS {} \;
	find $(distdir) -name "*.cc" -exec ./addauthors.sh ./AUTHORS {} \;
