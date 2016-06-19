#install headers
install(FILES
    capabilities.hh  datahandle.hh  entitypointer.hh  indexset.hh  intersectioniterator.hh
    entity.hh        geometry.hh    intersection.hh   iterator.hh
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/geometrygridpart)
