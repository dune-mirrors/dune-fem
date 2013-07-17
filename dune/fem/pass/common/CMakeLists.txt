set(passcommondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/pass/common)
set(passcommon_HEADERS  localfunctionselector.hh
                     localfunctiontuple.hh
                     local.hh
                     pass.hh
                     selector.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${passcommon_HEADERS} DESTINATION ${passcommondir})