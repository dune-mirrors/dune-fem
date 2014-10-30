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

set(TESTS lagrangeinterpolation_linear 
          lagrangeinterpolation_quadratic
          dgl2projection_adaptive
          dgl2projection_blockvector
          dgl2projection_vector
          dgl2projection_vector_int
          dgl2projection_attached
          dgl2projection_legendre
          dgl2projection_hierarchiclegendre
          dgl2projection_blockvectordiscretefunction
          dgl2projection_fv
          dgl2projection_idgridpart
          vtxprojection_adaptive
          intersectionindexset)

# non compiling tests
# todo: fix compiler errors
#set(TESTS "${TESTS} dgl2projection_filtered")


# copy data to build source to make tests work
configure_file(1dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/1dgrid.dgf COPYONLY)
configure_file(2dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/2dgrid.dgf COPYONLY)
configure_file(3dgrid.dgf ${CMAKE_CURRENT_BINARY_DIR}/3dgrid.dgf COPYONLY)
configure_file(parameter ${CMAKE_CURRENT_BINARY_DIR}/parameter COPYONLY)
configure_file(adaptation-4.sh ${CMAKE_CURRENT_BINARY_DIR}/adaptation_4.sh COPYONLY)

add_executable(lagrangeinterpolation_linear lagrangeinterpolation.cc)
dune_target_link_libraries(lagrangeinterpolation_linear "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeinterpolation_linear APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(lagrangeinterpolation_quadratic lagrangeinterpolation.cc)
dune_target_link_libraries(lagrangeinterpolation_quadratic "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET lagrangeinterpolation_quadratic APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=2;COUNT_FLOPS;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_adaptive dgl2projection.cc)
dune_target_link_libraries(dgl2projection_adaptive "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_adaptive APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

# only run this test for YaspGrid
add_executable(dgl2projection_legendre dgl2projection.cc)
dune_target_link_libraries(dgl2projection_legendre "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_legendre APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_LEGENDRESPACE;YASPGRID;GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_hierarchiclegendre dgl2projection.cc)
dune_target_link_libraries(dgl2projection_hierarchiclegendre "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_hierarchiclegendre APPEND PROPERTY 
  COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_HIERARCHICLEGENDRESPACE;YASPGRID;GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_filtered dgl2projection.cc)
dune_target_link_libraries(dgl2projection_filtered "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_filtered APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_FILTEREDGRID;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_fv dgl2projection.cc)
dune_target_link_libraries(dgl2projection_fv "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_fv APPEND PROPERTY 
             COMPILE_DEFINITIONS "COUNT_FLOPS;USE_FVSPACE;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_blockvector dgl2projection.cc)
dune_target_link_libraries(dgl2projection_blockvector "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_blockvector APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_BLOCKVECTORFUNCTION;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_vector dgl2projection.cc)
dune_target_link_libraries(dgl2projection_vector "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_vector  APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;USE_VECTORFUNCTION;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_vector_int dgl2projection.cc)
dune_target_link_libraries(dgl2projection_vector_int "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_vector_int APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;USE_VECTORFUNCTION;USE_DOFTYPE_INT;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_attached dgl2projection.cc)
dune_target_link_libraries(dgl2projection_attached "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_attached APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;USE_ATTACHEDFUNCTION;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_blockvectordiscretefunction dgl2projection.cc)
dune_target_link_libraries(dgl2projection_blockvectordiscretefunction "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_blockvectordiscretefunction APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;USE_BLOCKVECTORDISCRETEFUNCTION;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(dgl2projection_idgridpart dgl2projection.cc)
dune_target_link_libraries(dgl2projection_idgridpart "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET dgl2projection_idgridpart APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;USE_IDGRIDPART;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(vtxprojection_adaptive vtxprojection.cc)
dune_target_link_libraries(vtxprojection_adaptive "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET vtxprojection_adaptive APPEND PROPERTY 
             COMPILE_DEFINITIONS "POLORDER=1;COUNT_FLOPS;${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(intersectionindexset intersectionindexset.cc)
dune_target_link_libraries(intersectionindexset "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET intersectionindexset APPEND PROPERTY 
             COMPILE_DEFINITIONS "${GRIDTYPE};GRIDDIM=${GRIDDIM}" )

add_executable(adaptation adaptation.cc)
dune_target_link_libraries(adaptation "${DUNE_LIBS};${LOCAL_LIBS}")
set_property(TARGET adaptation APPEND PROPERTY 
  COMPILE_DEFINITIONS "GRIDDIM=3;WORLDDIM=3;${GRIDTYPE}" )

### example "make GRIDTYPE=..."-hack:
#set_gridtype_from_command_line(adaptation GRIDDIM ${GRIDDIM} WORLDDIM ${WORLDDIM} GRIDTYPE ${GRIDTYPE}) 
#set_gridtype_from_command_line(adaptation FORCE) 



set(ALL_TESTS ${TESTS} adaptation)

foreach(i ${ALL_TESTS})
  add_dune_zlib_flags(${i}) 
  ###create all known gridtype-griddim-worlddim combinations
  #dune_create_gridtype_executable(${i})
endforeach()

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${TESTS})
  add_test(${i} ${i})
endforeach(i ${TESTS})
add_test(adaptation_4.sh adaptation_4.sh)
# We do not want want to build the tests during make all
# but just build them on demand
add_directory_test_target(_test_target)
add_dependencies(${_test_target} ${ALL_TESTS} ${COMPILE_XFAIL_TESTS})

