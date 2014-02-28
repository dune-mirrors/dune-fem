set(matlabdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/matlab)

set(matlab_HEADERS  matlabhelper.hh
                 matlabstream.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${matlab_HEADERS} DESTINATION ${matlabdir})
