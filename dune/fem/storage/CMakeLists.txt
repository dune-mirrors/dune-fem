# don't follow the full GNU-standard; we definitely need automake >= 1.5
set(AUTOMAKE_OPTIONS  foreign 1.5)
set(storagedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/storage)
set(storage_HEADERS  arrayallocator.hh
                        array.hh
                                                                        array_inline.hh
                              blasvector.hh
                          envelope.hh
                              objectstack.hh
                          obstack.hh
                              rectarray.hh
                          referencecounter.hh
                              stlarray.hh
                          singletonlist.hh
                              subarray.hh
                          vector.hh
                              vector_inline.hh
                          ssematrixdynamic.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${storage_HEADERS} DESTINATION ${storagedir})