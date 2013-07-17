set(mapperdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/mapper)
set(mapper_HEADERS  codimensionmapper.hh  dofmapper.hh  genericdofmapper.hh
                 nonblockmapper.hh genericadaptivedofmapper.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${mapper_HEADERS} DESTINATION ${mapperdir})