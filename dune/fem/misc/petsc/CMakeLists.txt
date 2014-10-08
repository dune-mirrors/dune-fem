# don't follow the full GNU-standard; we definitely need automake >= 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)
set(petscdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/misc/petsc)
set(petsc_HEADERS petsccommon.hh
                  petscdofblock.hh
                  petscdofmappings.hh
                  petscghostarraybuilder.hh
                  petscslavedofprovider.hh
                  petscvector.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${petsc_HEADERS} DESTINATION ${petscdir})
