set(dofmapperdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/dofmapper)
set(dofmapper_HEADERS  code.hh  compile.hh  indexsetdofmapper.hh  exceptions.hh  localkey.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${dofmapper_HEADERS} DESTINATION ${dofmapperdir})