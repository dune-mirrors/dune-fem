set(projectiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator/projection)
set(projection_HEADERS  hdivprojection.hh l2projection.hh dgl2projection.hh vtxprojection.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${projection_HEADERS} DESTINATION ${projectiondir})