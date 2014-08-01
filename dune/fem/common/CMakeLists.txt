set(femcommondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/common)
set(femcommon_HEADERS  tupletypetraits.hh
                    tupleutility.hh
                    typeindexedtuple.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${femcommon_HEADERS} DESTINATION ${femcommondir})

