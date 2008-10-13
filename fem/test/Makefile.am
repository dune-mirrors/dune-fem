# IGNORE_HEADERS

EXTRA_DIST = 2dgrid.dgf \
	     3dgrid.dgf \
	     dfspace.hh \
	     dgl2projection.hh \
	     exactsolution.hh \
	     testgrid.hh

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
AM_CXXFLAGS = -O0 -g -Wall $(ALL_PKG_CXXFLAGS) $(MPI_CPPFLAGS)

TESTPROGS = lagrangeinterpolation_linear lagrangeinterpolation_quadratic \
            dgl2projection_adaptive      dgl2projection_blockvector \
	    dgl2projection_vector 	 dgl2projection_attached \
	    dgl2projection_fv \
	    dgl2projection_cmbpoint      dgl2projection_cmbvariable \
	    vtxprojection_adaptive \
	    mapper_dg                    mapper_fv \
	    mapper_cmbpoint              mapper_cmbvariable


TESTS = $(TESTPROGS)
check_PROGRAMS = $(TESTPROGS)

lagrangeinterpolation_linear_SOURCES = lagrangeinterpolation.cc
lagrangeinterpolation_linear_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

lagrangeinterpolation_quadratic_SOURCES = lagrangeinterpolation.cc
lagrangeinterpolation_quadratic_CPPFLAGS = -DPOLORDER=2 -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_adaptive_SOURCES = dgl2projection.cc
dgl2projection_adaptive_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_fv_SOURCES = dgl2projection.cc
dgl2projection_fv_CPPFLAGS = -DUSE_FVSPACE -DCOUNT_FLOPS $(AM_CPPFLAGS)

dgl2projection_cmbpoint_SOURCES = dgl2projection.cc
dgl2projection_cmbpoint_CPPFLAGS = -DPOLORDER=1 -DUSE_COMBINEDSPACE $(AM_CPPFLAGS)

dgl2projection_cmbvariable_SOURCES = dgl2projection.cc
dgl2projection_cmbvariable_CPPFLAGS = -DPOLORDER=1 -DUSE_COMBINEDSPACE -DUSE_VARIABLEBASE $(AM_CPPFLAGS)

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

mapper_fv_SOURCES = mapper.cc
mapper_fv_CPPFLAGS = -DUSE_FVSPACE -DCOUNT_FLOPS $(AM_CPPFLAGS)

mapper_cmbpoint_SOURCES = mapper.cc
mapper_cmbpoint_CPPFLAGS = -DPOLORDER=1 -DUSE_COMBINEDSPACE $(AM_CPPFLAGS)

mapper_cmbvariable_SOURCES = mapper.cc
mapper_cmbvariable_CPPFLAGS = -DPOLORDER=1 -DUSE_COMBINEDSPACE -DUSE_VARIABLEBASE $(AM_CPPFLAGS)

include $(top_srcdir)/am/global-rules

