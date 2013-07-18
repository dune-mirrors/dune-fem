set(SUBDIRS  caching test)
set(quadraturedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/quadrature)
set(quadrature_HEADERS  cachingpointlist.hh
                     cachingquadrature.hh
                     cornerpointset.hh
                     defaultquadratures.hh
                     dunequadratures.hh
                     elementpointlist.hh
                     elementpointlistbase.hh
                     elementquadrature.hh
                     femquadratures.hh
                     femquadratures_inline.hh
                     gausspoints.hh
                     gausspoints_implementation.hh
                     idprovider.hh
                     integrator.hh
                     quadrature.hh
                     quadratureimp.hh
                     quadratureimp_inline.hh
                     pyramidpoints.hh
                     quadprovider.hh
                     simplexpoints.hh
                     intersectionquadrature.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${quadrature_HEADERS} DESTINATION ${quadraturedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})

dune_add_library(quadrature SOURCES 
    caching/twistutility.cc 
    caching/topology.cc
  )
