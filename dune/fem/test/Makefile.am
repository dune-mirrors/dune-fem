# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
GRIDDIM=2

EXTRA_DIST = 2dgrid.dgf \
             3dgrid.dgf \
             dfspace.hh \
             exactsolution.hh \
             testgrid.hh

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILIBS) $(LOCAL_LIBS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
AM_CXXFLAGS = -O0 -g -Wall $(ALL_PKG_CXXFLAGS) $(DUNEMPICPPFLAGS)

TESTPROGS = lagrangeinterpolation_linear lagrangeinterpolation_quadratic \
            dgl2projection_adaptive      dgl2projection_filtered \
            dgl2projection_blockvector   dgl2projection_vector \
            dgl2projection_attached \
            dgl2projection_fv      \
            vtxprojection_adaptive \
            mapper_dg                    mapper_fv    mapper_lagrange


TESTS = $(TESTPROGS)
check_PROGRAMS = $(TESTPROGS)

lagrangeinterpolation_linear_SOURCES = lagrangeinterpolation.cc
lagrangeinterpolation_linear_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

lagrangeinterpolation_quadratic_SOURCES = lagrangeinterpolation.cc
lagrangeinterpolation_quadratic_CPPFLAGS = -DPOLORDER=2 -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_adaptive_SOURCES = dgl2projection.cc
dgl2projection_adaptive_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_filtered_SOURCES = dgl2projection.cc
dgl2projection_filtered_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_FILTEREDGRID 

dgl2projection_fv_SOURCES = dgl2projection.cc
dgl2projection_fv_CPPFLAGS = -DUSE_FVSPACE -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_blockvector_SOURCES = dgl2projection.cc
dgl2projection_blockvector_CPPFLAGS = -DPOLORDER=1 -DUSE_BLOCKVECTORFUNCTION $(AM_CPPFLAGS)

dgl2projection_vector_SOURCES = dgl2projection.cc
dgl2projection_vector_CPPFLAGS = -DPOLORDER=1 -DUSE_VECTORFUNCTION $(AM_CPPFLAGS)

dgl2projection_attached_SOURCES = dgl2projection.cc
dgl2projection_attached_CPPFLAGS = -DPOLORDER=1 -DUSE_ATTACHEDFUNCTION $(AM_CPPFLAGS)

vtxprojection_adaptive_SOURCES = vtxprojection.cc
vtxprojection_adaptive_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

mapper_dg_SOURCES = mapper.cc
mapper_dg_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

mapper_lagrange_SOURCES = mapper.cc
mapper_lagrange_CPPFLAGS = -DUSE_LAGRANGESPACE -DPOLORDER=2 -DCOUNT_FLOPS $(AM_CPPFLAGS)

mapper_fv_SOURCES = mapper.cc
mapper_fv_CPPFLAGS = -DUSE_FVSPACE -DCOUNT_FLOPS $(AM_CPPFLAGS)

include $(top_srcdir)/am/global-rules

