set(discontinuousgalerkindir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/discontinuousgalerkin)
set(discontinuousgalerkin_HEADERS  capabilities.hh
                                declaration.hh
                                default.hh
                                lagrange.hh
                                legendre.hh
                                legendrepolynomials.hh
                                localdgmassmatrix.hh
                                localinterpolation.hh
                                localrestrictprolong.hh
                                orthonormalbase_mod.hh
                                shapefunctionset.hh
                                space.hh)
set(EXTRA_DIST  orthonormalbase_1d.cc orthonormalbase_2d.cc orthonormalbase_3d.cc)
set(noinst_LTLIBRARIES  libdiscontinuousgalerkin.la)

set(libdiscontinuousgalerkin_la_LIBADD )
dune_add_library(discontinuousgalerkin SOURCES legendrepolynomials.cc
                                      orthonormalbase_1d.cc
                                      orthonormalbase_2d.cc
                                      orthonormalbase_3d.cc)
target_link_libraries(discontinuousgalerkin ${DUNE_LIBS})
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${discontinuousgalerkin_HEADERS} DESTINATION ${discontinuousgalerkindir})