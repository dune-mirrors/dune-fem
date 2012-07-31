parallel_laplacedir = $(includedir)/dune/fem/petsc/experiments/parallel_laplace

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILIBS) $(LOCAL_LIBS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

if PETSC
CHECK_PROGS = parallel_laplace
else
CHECK_PROGS = 
endif

check_PROGRAMS = $(CHECK_PROGS)

CLEANFILES = *.tex data/* rhs.out.gnuplot* solution.out.gnuplot* error.out.gnuplot*

parallel_laplace_SOURCES = parallel_laplace.cc
parallel_laplace_CPPFLAGS = $(AM_CPPFLAGS)

include $(top_srcdir)/am/global-rules
