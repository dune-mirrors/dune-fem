set(SUBDIRS  test)
set(commondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator/common)
set(common_HEADERS  automaticdifferenceoperator.hh  differentiableoperator.hh
                 discreteoperator.hh  localoperatordef.inc  operator.hh
                 combine.inc  localoperator.hh  mapping.hh  objpointer.hh
                 vectorspace.hh  spaceoperatorif.hh  localmatrix.hh
                 localmatrixwrapper.hh  temporarylocalmatrix.hh
                                                                 petsclinearoperator.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${common_HEADERS} DESTINATION ${commondir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
