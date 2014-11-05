# $Id: Makefile.am 6360 2011-02-28 11:41:46Z sander $
set(SUBDIRS  pages)
# setting like in dune-web
set(BASEDIR ../..)
set(CURDIR doc/doxygen)
# add some more stuff to install and tarball
set(DOCFILES  ${PAGES})
set(EXTRA_DIST  ${PAGES})
# include not needed for CMake
# include $(top_srcdir)/am/doxygen
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
dune_add_subdirs(${SUBDIRS})
# Install documentation
# We assume that ${DOCFILES} and ${DOCFILES_EXTRA} are targets
install(TARGETS ${DOCFILES} ${DOCFILES_EXTRA} DESTINATION ${docdir})
# Create Doxyfile.in and Doxyfile, and doxygen documentation
add_doxygen_target()
