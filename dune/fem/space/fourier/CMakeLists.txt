set(fourierdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/fourier)
set(fourier_HEADERS  capabilities.hh
                  declaration.hh
                  dofmapper.hh
                  functionset.hh
                  space.hh)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${fourier_HEADERS} DESTINATION ${fourierdir})
