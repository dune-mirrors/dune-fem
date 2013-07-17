set(fourierdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/fourier)
set(fourier_HEADERS  capabilities.hh
                  declaration.hh
                  dofmapper.hh
                  functionset.hh
                  space.hh)
set(SUBDIRS )
set(EXTRA_DIST )
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${fourier_HEADERS} DESTINATION ${fourierdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
