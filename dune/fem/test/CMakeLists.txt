# IGNORE_HEADERS
# set default grid type for make all
set(GRIDTYPE YASPGRID)
#GRIDTYPE=ALUGRID_CONFORM
set(GRIDDIM 2)
#GRIDDIM=3
set(EXTRA_DIST  2dgrid.dgf
             3dgrid.dgf
                                                 adaptation-4.sh
             dfspace.hh
             exactsolution.hh
             weightfunction.hh
             testgrid.hh)
# linker flags
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${DUNEMPILIBS} ${LOCAL_LIBS})
# preprocessor flags
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS})
set(TESTPROGS  lagrangeinterpolation_linear lagrangeinterpolation_quadratic
            dgl2projection_adaptive      dgl2projection_filtered
            dgl2projection_blockvector   dgl2projection_vector dgl2projection_vector_int
            dgl2projection_attached      dgl2projection_legendre
            dgl2projection_blockvectordiscretefunction
            dgl2projection_fv
            dgl2projection_idgridpart
            vtxprojection_adaptive
            intersectionindexset)
set(TESTS  ${TESTPROGS} adaptation-4.sh)
set(check_PROGRAMS  ${TESTPROGS} adaptation)
add_executable(lagrangeinterpolation_linear lagrangeinterpolation.cc)
target_link_libraries(lagrangeinterpolation_linear ${DUNE_LIBS})
add_executable(lagrangeinterpolation_quadratic lagrangeinterpolation.cc)
target_link_libraries(lagrangeinterpolation_quadratic ${DUNE_LIBS})
add_executable(dgl2projection_adaptive dgl2projection.cc)
target_link_libraries(dgl2projection_adaptive ${DUNE_LIBS})
# only run this test for YaspGrid
add_executable(dgl2projection_legendre dgl2projection.cc)
target_link_libraries(dgl2projection_legendre ${DUNE_LIBS})
add_executable(dgl2projection_filtered dgl2projection.cc)
target_link_libraries(dgl2projection_filtered ${DUNE_LIBS})
add_executable(dgl2projection_fv dgl2projection.cc)
target_link_libraries(dgl2projection_fv ${DUNE_LIBS})
add_executable(dgl2projection_blockvector dgl2projection.cc)
target_link_libraries(dgl2projection_blockvector ${DUNE_LIBS})
add_executable(dgl2projection_vector dgl2projection.cc)
target_link_libraries(dgl2projection_vector ${DUNE_LIBS})
add_executable(dgl2projection_vector_int dgl2projection.cc)
target_link_libraries(dgl2projection_vector_int ${DUNE_LIBS})
add_executable(dgl2projection_attached dgl2projection.cc)
target_link_libraries(dgl2projection_attached ${DUNE_LIBS})
add_executable(dgl2projection_blockvectordiscretefunction dgl2projection.cc)
target_link_libraries(dgl2projection_blockvectordiscretefunction ${DUNE_LIBS})
add_executable(dgl2projection_idgridpart dgl2projection.cc)
target_link_libraries(dgl2projection_idgridpart ${DUNE_LIBS})
add_executable(vtxprojection_adaptive vtxprojection.cc)
target_link_libraries(vtxprojection_adaptive ${DUNE_LIBS})
add_executable(intersectionindexset intersectionindexset.cc)
target_link_libraries(intersectionindexset ${DUNE_LIBS})
add_executable(adaptation adaptation.cc)
target_link_libraries(adaptation ${DUNE_LIBS})
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})