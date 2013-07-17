set(idgridpartdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/idgridpart)
set(idgridpart_HEADERS  capabilities.hh  datahandle.hh  entity.hh
                     entitypointer.hh  geometry.hh  indexset.hh
                     intersection.hh  intersectioniterator.hh
                     iterator.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${idgridpart_HEADERS} DESTINATION ${idgridpartdir})