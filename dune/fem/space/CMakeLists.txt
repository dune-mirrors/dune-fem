set(spacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space)
set(space_HEADERS  combineddiscretefunctionspace.hh
                combinedspace.hh
                discontinuousgalerkin.hh
                finitevolume.hh
                fourier.hh
                lagrange.hh
                padaptivespace.hh
                rannacherturek.hh)
set(SUBDIRS  basisfunctionset
             combinedspace 
             common
             discontinuousgalerkin 
             dofmapper
             finitevolume
             fourier
             lagrange
             mapper
             padaptivespace 
             rannacherturek
             shapefunctionset
             test)

install(FILES ${space_HEADERS} DESTINATION ${spacedir})

dune_add_subdirs(${SUBDIRS})

dune_add_library(space OBJECT 
    shapefunctionset/legendrepolynomials.cc
    shapefunctionset/orthonormal/orthonormalbase_1d.cc
    shapefunctionset/orthonormal/orthonormalbase_2d.cc
    shapefunctionset/orthonormal/orthonormalbase_3d.cc
    lagrange/shapefunctionset.cc
  )
