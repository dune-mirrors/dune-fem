set(attachedfunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/attachedfunction)
set(attachedfunction_HEADERS  container.hh
                           function.hh
                           slotiterator.hh
                           storage.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${attachedfunction_HEADERS} DESTINATION ${attachedfunctiondir})