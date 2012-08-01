l2projectiondir = $(includedir)/dune/fem/petsc/tests/l2projection

#GRIDTYPE = SPGRID
GRIDTYPE = YASPGRID
GRIDDIM = 2
POLORDER = 1

EXTRA_DIST = parameter unitcube1.dgf unitcube2.dgf unitcube3.dgf base.hh massoperator.hh \
						 mpiexec-petscmatrix.sh mpiexec.sh 
	

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILIBS) $(LOCAL_LIBS) $(MPI_LDFLAGS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

if PETSC
CHECK_PROGS = l2projection l2projection_petscmatrix
else
CHECK_PROGS = 
endif

check_PROGRAMS = $(CHECK_PROGS)

CLEANFILES = *.tex data/*

l2projection_SOURCES = l2projection.cc
l2projection_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER)

l2projection_petscmatrix_SOURCES = l2projection.cc
l2projection_petscmatrix_CPPFLAGS = $(AM_CPPFLAGS) -DPOLORDER=$(POLORDER) -DPETSCLINEAROPERATOR=1

include $(top_srcdir)/am/global-rules
