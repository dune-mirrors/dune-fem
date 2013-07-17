# $Id: $
#the dune-fem library
set(lib_LTLIBRARIES  libdunefem.la)
# this is just a renaming of libfem to prevent name clashes
add_executable(libdunefem_la )
set(libdunefem_la_LIBADD  ${top_srcdir}/dune/fem/libfem.la)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
