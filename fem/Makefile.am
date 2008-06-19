# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# do everything in the subdirs
SUBDIRS = pass space function gridpart quadrature operator solver storage io misc test

noinst_LTLIBRARIES = libfem.la

# construct only from sub-libraries
libfem_la_SOURCES =
libfem_la_LIBADD = space/libspace.la solver/libsolver.la 

# install all files to dune/fem 
femdir = $(includedir)/dune/fem
fem_HEADERS = version.hh

include $(top_srcdir)/am/global-rules
