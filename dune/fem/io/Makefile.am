# install header files
iodir =  $(includedir)/dune/fem/io
io_HEADERS = io.hh  parameter.hh

# do everything in the subdirs
SUBDIRS = file streams visual matlab

# build io library
noinst_LTLIBRARIES = libio.la

libio_la_SOURCES = io.cc

# include common rules
include $(top_srcdir)/am/global-rules
