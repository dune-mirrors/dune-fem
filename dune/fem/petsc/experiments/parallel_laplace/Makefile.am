parallel_laplacedir = $(includedir)/dune/fem/petsc/experiments/parallel_laplace
parallel_laplace_HEADERS= function.hh l2error.hh laplaceassembler.hh printer.hh \
												  scheme.hh solver.hh

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

EXTRA_DIST = mpiexec.sh 

include $(top_srcdir)/am/global-rules
