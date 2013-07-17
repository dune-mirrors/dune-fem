set(localfunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/localfunction)
set(localfunction_HEADERS  average.hh
                        default.hh
                        functor.hh
                        generic.hh
                        generic_inline.hh
                        localfunction.hh
                        localfunctionsetadapter.hh
                        localfunctionset.hh
                        standard.hh
                        temporary.hh
                        temporarylocalfunction.hh
                        wrapper.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${localfunction_HEADERS} DESTINATION ${localfunctiondir})