# $Id$

# we need the module file to be able to build via dunecontrol
EXTRA_DIST=dune.module

# don't follow the full GNU-standard; we definitely need automake >= 1.5
#AUTOMAKE_OPTIONS = foreign 1.5

# do everything in the subdirs
DIST_SUBDIRS = fem lib doc
if BUILD_DOCS
SUBDIRS = fem lib doc
else
SUBDIRS = fem lib
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = dune-fem.pc

DISTCHECK_CONFIGURE_FLAGS = --with-dune="$(DUNEROOT)" CXX="$(CXX)" CC="$(CC)"
#DISTCHECK_CONFIGURE_FLAGS = --with-dune-common="$(DUNE_COMMON_ROOT)" --with-dune-grid="$(DUNE_GRID_ROOT)" CXX="$(CXX)" CC="$(CC)" DOXYGEN="true"

include $(top_srcdir)/am/global-rules

SOURCEFILE_MASK = -name \*.h -o -name \*.hh -o -name \*.c -o -name \*.cc

# add library name, authors list, and license to progam files
dist-hook:
	find $(distdir) \( $(SOURCEFILE_MASK) \) -exec scripts/addauthors.sh ./AUTHORS {} \;
