set(padaptivespacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/padaptivespace)
set(padaptivespace_HEADERS  adaptmanager.hh
                         capabilities.hh
                         declaration.hh
                         discontinuousgalerkin.hh
                         generic.hh
                         lagrange.hh
                         mapper.hh
                         padaptivespace.hh
                         restrictprolong.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${padaptivespace_HEADERS} DESTINATION ${padaptivespacedir})
