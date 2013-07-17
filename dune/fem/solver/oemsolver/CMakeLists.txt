set(oemsolverdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/solver/oemsolver)
set(oemsolver_HEADERS  bicgsq.h  bicgstab.h  cblas.h  cghs.h  gmres.h
oemsolver.hh  preconditioning.hh  tmpmem.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${oemsolver_HEADERS} DESTINATION ${oemsolverdir})