set(matlabdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/matlab)
set(SUBDIRS  test)
set(DIST_SUBDIRS  m test)
set(matlab_HEADERS  matlabhelper.hh
                 matlabstream.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${matlab_HEADERS} DESTINATION ${matlabdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
