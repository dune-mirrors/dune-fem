set(combinedspacedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/combinedspace)
set(combinedspace_HEADERS  combinedadaptmanager.hh
                        combineddofstorage.hh
                        combinedspace.hh
                        mapper.hh
                        combineddiscretefunctionspace.hh
                        combinedmapper.hh
                        combinedbasisfunctionset.hh
                        adaptmanager.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${combinedspace_HEADERS} DESTINATION ${combinedspacedir})