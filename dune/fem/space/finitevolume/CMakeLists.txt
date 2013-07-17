set(finitevolumespacedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/finitevolume)
set(finitevolumespace_HEADERS  declaration.hh  space.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${finitevolumespace_HEADERS} DESTINATION ${finitevolumespacedir})