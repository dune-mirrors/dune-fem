set(shapefunctionsetdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/shapefunctionset)
set(shapefunctionset_HEADERS  caching.hh
                           lagrange.hh
                           legendre.hh
                           localfunctions.hh
                           orthonormal.hh
                           proxy.hh
                           selectcaching.hh
                           shapefunctionset.hh
                           simple.hh
                           tensorproduct.hh
                           vectorial.hh
                           wrapper.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${shapefunctionset_HEADERS} DESTINATION ${shapefunctionsetdir})