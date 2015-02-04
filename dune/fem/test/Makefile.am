# IGNORE_HEADERS

# set default grid type for make all
GRIDTYPE=YASPGRID
#GRIDTYPE=ALUGRID_CONFORM

GRIDDIM=2
#GRIDDIM=3

EXTRA_DIST = 2dgrid.dgf \
             3dgrid.dgf \
						 adaptation-4.sh \
             dfspace.hh \
             exactsolution.hh \
						 massoperator.hh \
             weightfunction.hh \
             testgrid.hh

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILIBS) $(LOCAL_LIBS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)

TESTPROGS = combinedspaceadaptation \
						dgl2projection_filtered \
						dgl2projection_hierarchiclegendre\
						dgl2projection_legendre \
						dgl2projection_vector	\
						dgl2projection_vector_int \
						fmatrixconverter \
						l2projection_adaptive \
						l2projection_istl \
						l2projection_petsc \
						lagrangeinterpolation_linear \
						lagrangeinterpolation_quadratic \
						dgl2projection_adaptive	\
						dgl2projection_blockvector \
						dgl2projection_blockvectordiscretefunction \
						dgl2projection_fv      \
						dgl2projection_idgridpart \
						dgl2projection_petscdiscretefunction \
						intersectionindexset \
						vtxprojection_adaptive

TESTS = $(TESTPROGS) adaptation-4.sh
check_PROGRAMS = $(TESTPROGS) adaptation


adaptation_CPPFLAGS = $(AM_CPPFLAGS) -UGRIDDIM -DGRIDDIM=3 -UWORLDDIM -DWORLDDIM=3
adaptation_SOURCES = adaptation.cc

combinedspaceadaptation_CPPFLAGS = $(AM_CPPFLAGS) -DUSECOMBINEDSPACE -UGRIDDIM -DGRIDDIM=3 -UWORLDDIM -DWORLDDIM=3
combinedspaceadaptation_SOURCES = adaptation.cc

dgl2projection_adaptive_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)
dgl2projection_adaptive_SOURCES = dgl2projection.cc

dgl2projection_blockvector_CPPFLAGS = -DPOLORDER=1 -DUSE_BLOCKVECTORFUNCTION $(AM_CPPFLAGS)
dgl2projection_blockvector_SOURCES = dgl2projection.cc

dgl2projection_blockvectordiscretefunction_CPPFLAGS = -DPOLORDER=1 -DUSE_BLOCKVECTORDISCRETEFUNCTION $(AM_CPPFLAGS)
dgl2projection_blockvectordiscretefunction_SOURCES = dgl2projection.cc

dgl2projection_filtered_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_FILTEREDGRID 
dgl2projection_filtered_SOURCES = dgl2projection.cc

dgl2projection_fv_CPPFLAGS = -DUSE_FVSPACE -DCOUNT_FLOPS $(AM_CPPFLAGS)
dgl2projection_fv_SOURCES = dgl2projection.cc

# only run this test for YaspGrid
dgl2projection_hierarchiclegendre_CPPFLAGS = -DPOLORDER=1 -DUSE_HIERARCHICLEGENDRESPACE -DCOUNT_FLOPS $(AM_CPPFLAGS) -U$(GRIDTYPE) -DYASPGRID
dgl2projection_hierarchiclegendre_SOURCES = dgl2projection.cc

dgl2projection_idgridpart_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS) -DUSE_IDGRIDPART 
dgl2projection_idgridpart_SOURCES = dgl2projection.cc

# only run this test for YaspGrid
dgl2projection_legendre_CPPFLAGS = -DPOLORDER=1 -DUSE_LEGENDRESPACE -DCOUNT_FLOPS $(AM_CPPFLAGS) -U$(GRIDTYPE) -DYASPGRID
dgl2projection_legendre_SOURCES = dgl2projection.cc

dgl2projection_petscdiscretefunction_CPPFLAGS = -DPOLORDER=1 -DUSE_PETSCDISCRETEFUNCTION $(AM_CPPFLAGS)
dgl2projection_petscdiscretefunction_SOURCES = dgl2projection.cc

dgl2projection_vector_CPPFLAGS = -DPOLORDER=1 -DUSE_VECTORFUNCTION $(AM_CPPFLAGS)
dgl2projection_vector_SOURCES = dgl2projection.cc

dgl2projection_vector_int_CPPFLAGS = -DPOLORDER=1 -DUSE_VECTORFUNCTION -DUSE_DOFTYPE_INT $(AM_CPPFLAGS)
dgl2projection_vector_int_SOURCES = dgl2projection.cc

fmatrixconverter_CPPFLAGS = -DCOUNT_FLOPS $(AM_CPPFLAGS)
fmatrixconverter_SOURCES = fmatrixconverter.cc

intersectionindexset_CPPFLAGS = $(AM_CPPFLAGS)
intersectionindexset_SOURCES = intersectionindexset.cc

l2projection_adaptive_CPPFLAGS = $(AM_CPPFLAGS) -DCOUNT_FLOPS -DPOLORDER=1
l2projection_adaptive_SOURCES = l2projection.cc

l2projection_istl_CPPFLAGS = $(AM_CPPFLAGS) -DCOUNT_FLOPS -DPOLORDER=1 -DUSE_BLOCKVECTORFUNCTION
l2projection_istl_SOURCES = l2projection.cc

l2projection_petsc_CPPFLAGS = $(AM_CPPFLAGS) -DCOUNT_FLOPS -DPOLORDER=1 -DUSE_PETSCDISCRETEFUNCTION
l2projection_petsc_SOURCES = l2projection.cc

lagrangeinterpolation_linear_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)
lagrangeinterpolation_linear_SOURCES = lagrangeinterpolation.cc

lagrangeinterpolation_quadratic_CPPFLAGS = -DPOLORDER=2 -DCOUNT_FLOPS $(AM_CPPFLAGS)
lagrangeinterpolation_quadratic_SOURCES = lagrangeinterpolation.cc

vtxprojection_adaptive_CPPFLAGS = -DPOLORDER=1 -DCOUNT_FLOPS $(AM_CPPFLAGS)
vtxprojection_adaptive_SOURCES = vtxprojection.cc

include $(top_srcdir)/am/global-rules
