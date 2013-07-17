set(femcommondir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/common)
set(femcommon_HEADERS  tupletypetraits.hh
                    tupleutility.hh
                    typeindexedtuple.hh)
set(SUBDIRS )
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${femcommon_HEADERS} DESTINATION ${femcommondir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
