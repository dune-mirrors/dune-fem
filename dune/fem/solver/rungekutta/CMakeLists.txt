set(rungekuttadir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/solver/rungekutta)


set(rungekutta_HEADERS basicimplicit.hh
                       butchertable.hh
                       explicit.hh
                       implicit.hh
                       semiimplicit.hh
                       timestepcontrol.hh)

install(FILES ${oemsolver_HEADERS} DESTINATION ${oemsolverdir})
