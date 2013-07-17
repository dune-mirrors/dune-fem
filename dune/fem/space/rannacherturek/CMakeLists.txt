set(rannacherturekspacedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/rannacherturek)
set(rannacherturekspace_HEADERS  capabilities.hh
                              declaration.hh
                              dofmappercode.hh
                              localinterpolation.hh
                              space.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${rannacherturekspace_HEADERS} DESTINATION ${rannacherturekspacedir})