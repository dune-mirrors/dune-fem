# IGNORE_HEADERS

GRIDTYPE = YASPGRID
GRIDDIM = 2

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(MPI_LDFLAGS) ../libfem.la
# preprocessor flags
AM_CPPFLAGS = -Wall $(ALL_PKG_CPPFLAGS)

TESTPROGS = lagrangeinterpolation_linear lagrangeinterpolation_quadratic \
            dgl2projection_adaptive      dgl2projection_blockvector \
						dgl2projection_cmbpoint      dgl2projection_cmbvariable


TESTS = $(TESTPROGS)
check_PROGRAMS = $(TESTPROGS)

lagrangeinterpolation_linear_SOURCES = lagrangeinterpolation.cc
lagrangeinterpolation_linear_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

lagrangeinterpolation_quadratic_SOURCES = lagrangeinterpolation.cc
lagrangeinterpolation_quadratic_CPPFLAGS = -DPOLORDER=2 -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_adaptive_SOURCES = dgl2projection.cc
dgl2projection_adaptive_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_cmbpoint_SOURCES = dgl2projection.cc
dgl2projection_cmbpoint_CPPFLAGS = -DPOLORDER=1 -DUSE_COMBINEDSPACE $(AM_CPPFLAGS)

dgl2projection_cmbvariable_SOURCES = dgl2projection.cc
dgl2projection_cmbvariable_CPPFLAGS = -DPOLORDER=1 -DUSE_COMBINEDSPACE -DUSE_VARIABLEBASE $(AM_CPPFLAGS)

dgl2projection_blockvector_SOURCES = dgl2projection.cc
dgl2projection_blockvector_CPPFLAGS = -DPOLORDER=1 -DUSE_BLOCKVECTORFUNCTION $(AM_CPPFLAGS)
