set(commondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/common)
set(common_HEADERS  discretefunction.hh
                 discretefunction_inline.hh
                 dofblock.hh
                 dofiterator.hh
                 function.hh
                 functionset.hh
                 gridfunctionadapter.hh
                 localfunctionadapter.hh
                 scalarproducts.hh
                                                                 timedependentfunction.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${common_HEADERS} DESTINATION ${commondir})