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
set(EXTRA_DIST  orthonormalbase_1d.cc orthonormalbase_2d.cc orthonormalbase_3d.cc)

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${discontinuousgalerkin_HEADERS} DESTINATION ${discontinuousgalerkindir})
