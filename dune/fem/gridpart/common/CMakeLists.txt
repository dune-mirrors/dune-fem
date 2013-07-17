set(commondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/gridpart/common)
set(common_HEADERS  capabilities.hh  deaditerator.hh  gridpart.hh  gridpartview.hh defaultgridpartentity.hh metatwistutility.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${common_HEADERS} DESTINATION ${commondir})