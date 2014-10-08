set(localfunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/localfunction)
set(localfunction_HEADERS average.hh
                          const.hh
                          functor.hh
                          localfunction.hh
                          localfunctionsetadapter.hh
                          localfunctionset.hh
                          mutable.hh
                          temporary.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${localfunction_HEADERS} DESTINATION ${localfunctiondir})
