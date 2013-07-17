# only build this stuff if GRAPE was found
# einfach alle Pakete verwenden, die überhaupt gefunden wurden
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(datadispdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/visual/grape/datadisp)
set(datadisp_HEADERS  dataconvert.cc
                   datadisp.cc
                   errordisplay.hh
                   grcommon.hh
                   printhelp.cc
                   programtemplate.cc
                   readiotupledata.cc
                   readioparams.cc
                   readtupledata.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${datadisp_HEADERS} DESTINATION ${datadispdir})