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
