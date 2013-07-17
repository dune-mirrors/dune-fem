set(localdgpassdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/pass/localdg)
set(localdgpass_HEADERS  discretemodel.hh
                      modelcaller.hh
                      pass.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${localdgpass_HEADERS} DESTINATION ${localdgpassdir})