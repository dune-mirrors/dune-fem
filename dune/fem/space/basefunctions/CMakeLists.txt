set(basefunctionsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/basefunctions)
set(basefunctions_HEADERS  basefunctionfactory.hh  basefunctioninterface.hh
  basefunctionproxy.hh  basefunctionsetinterface.hh  basefunctionsets.hh
  basefunctionstorage.hh  basefunctionstorage.cc  genericbasefunctionsets.hh
  evaluatecaller.hh codegen.hh vectorialbasefunctionset.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${basefunctions_HEADERS} DESTINATION ${basefunctionsdir})