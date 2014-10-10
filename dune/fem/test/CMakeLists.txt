# IGNORE_HEADERS
# set default grid type for make all
set( GRIDTYPE YASPGRID )
#set( GRIDTYPE ALUGRID_CONFORM )
set(GRIDDIM 2 )
#set( GRIDDIM 3 )
set(EXTRA_DIST  2dgrid.dgf
                3dgrid.dgf
                adaptation-4.sh
                dfspace.hh
                exactsolution.hh
                weightfunction.hh
                testgrid.hh)

add_definitions( "-D${GRIDTYPE}" )
add_definitions( "-DGRIDDIM=${GRIDDIM}" ) 
           
set(TESTS adaptation_4.sh
          lagrangeinterpolation_linear 
          lagrangeinterpolation_quadratic
          dgl2projection_adaptive
          dgl2projection_filtered
          dgl2projection_blockvector
          dgl2projection_vector
          dgl2projection_vector_int
          dgl2projection_attached
          dgl2projection_legendre
          dgl2projection_blockvectordiscretefunction
          dgl2projection_fv
          dgl2projection_idgridpart
          vtxprojection_adaptive
          intersectionindexset)

add_executable(lagrangeinterpolation_linear lagrangeinterpolation.cc)
dune_target_link_libraries(lagrangeinterpolation_linear "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(lagrangeinterpolation_quadratic lagrangeinterpolation.cc)
dune_target_link_libraries(lagrangeinterpolation_quadratic "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_adaptive dgl2projection.cc)
dune_target_link_libraries(dgl2projection_adaptive "${DUNE_LIBS};${LOCAL_LIBS}")
# only run this test for YaspGrid
add_executable(dgl2projection_legendre dgl2projection.cc)
dune_target_link_libraries(dgl2projection_legendre "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_filtered dgl2projection.cc)
dune_target_link_libraries(dgl2projection_filtered "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_fv dgl2projection.cc)
dune_target_link_libraries(dgl2projection_fv "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_blockvector dgl2projection.cc)
dune_target_link_libraries(dgl2projection_blockvector "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_vector dgl2projection.cc)
dune_target_link_libraries(dgl2projection_vector "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_vector_int dgl2projection.cc)
dune_target_link_libraries(dgl2projection_vector_int "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_attached dgl2projection.cc)
dune_target_link_libraries(dgl2projection_attached "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_blockvectordiscretefunction dgl2projection.cc)
dune_target_link_libraries(dgl2projection_blockvectordiscretefunction "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(dgl2projection_idgridpart dgl2projection.cc)
dune_target_link_libraries(dgl2projection_idgridpart "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(vtxprojection_adaptive vtxprojection.cc)
dune_target_link_libraries(vtxprojection_adaptive "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(intersectionindexset intersectionindexset.cc)
dune_target_link_libraries(intersectionindexset "${DUNE_LIBS};${LOCAL_LIBS}")
add_executable(adaptation adaptation.cc)
dune_target_link_libraries(adaptation "${DUNE_LIBS};${LOCAL_LIBS}")

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${TESTS} ${COMPILE_XFAIL_TESTS})

