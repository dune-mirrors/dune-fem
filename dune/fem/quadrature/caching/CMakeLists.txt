set(SUBDIRS  test)
set(cachingdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/quadrature/caching)
set(caching_HEADERS  cacheprovider.cc  pointmapper.hh    pointprovider.hh  topology.cc  twistprovider.cc  twistutility.cc
                                                                        cacheprovider.hh  pointprovider.cc  registry.hh       topology.hh  twistprovider.hh  twistutility.hh)

install(FILES ${caching_HEADERS} DESTINATION ${cachingdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
