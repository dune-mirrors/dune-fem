set(discontinuousgalerkindir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/discontinuousgalerkin)
set(discontinuousgalerkin_HEADERS basisfunctionsets.hh
                                  declaration.hh
                                  generic.hh
                                  hierarchiclegendre.hh
                                  hierarchiclegendremap.hh
                                  interpolation.hh
                                  lagrange.hh
                                  legendre.hh
                                  localdgmassmatrix.hh
                                  localinterpolation.hh
                                  localrestrictprolong.hh
                                  shapefunctionsets.hh
                                  space.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${discontinuousgalerkin_HEADERS} DESTINATION ${discontinuousgalerkindir})
