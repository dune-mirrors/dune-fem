petsc_matrixdir = $(includedir)/dune/fem/petsc/experiments/petsc_matrix
petsc_matrix_HEADERS = function.hh l2error.hh laplaceassembler.hh printer.hh \
											 scheme.hh solver.hh

#EXTRA_DIST = parameter unitcube1.dgf unitcube2.dgf unitcube3.dgf

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILIBS) $(LOCAL_LIBS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

if PETSC
CHECK_PROGS = matrix_test
else
CHECK_PROGS = 
endif

check_PROGRAMS = $(CHECK_PROGS)

CLEANFILES = *.tex data/*

matrix_test_SOURCES = matrix_test.cc
matrix_test_CPPFLAGS = $(AM_CPPFLAGS)


include $(top_srcdir)/am/global-rules
