# do everything in the subdirs
SUBDIRS = pass space function gridpart quadrature operator solver storage io misc test

noinst_LTLIBRARIES = libfem.la

# construct only from sub-libraries
libfem_la_SOURCES =
libfem_la_LIBADD = space/libspace.la solver/libsolver.la quadrature/libquadrature.la

# install all files to dune/fem 
femdir = $(includedir)/dune/fem
fem_HEADERS = version.hh

include $(top_srcdir)/am/global-rules
