set(spacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space)
set(space_HEADERS  combineddiscretefunctionspace.hh
                combinedspace.hh
                dgspace.hh
                discontinuousgalerkin.hh
                finitevolume.hh
                fourier.hh
                fvspace.hh
                lagrange.hh
                lagrangespace.hh
                padaptivespace.hh
                rannacherturek.hh
                reducedbasisspace.hh)
set(SUBDIRS  basefunctions basisfunctionset
          combinedspace common
          discontinuousgalerkin dofmapper
          finitevolume fourier lagrange mapper
          padaptivespace rannacherturek
          reducedbasisspace shapefunctionset
          test)

install(FILES ${space_HEADERS} DESTINATION ${spacedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

dune_add_library(space SOURCES 
    discontinuousgalerkin/legendrepolynomials.cc
    discontinuousgalerkin/orthonormalbase_1d.cc
    discontinuousgalerkin/orthonormalbase_2d.cc
    discontinuousgalerkin/orthonormalbase_3d.cc
    lagrange/shapefunctionset.cc
  )