set(SUBDIRS  common localdg test)
set(passdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/pass)
set(pass_HEADERS  applylocaloperator.hh
               createpass.hh
               dgdiscretemodel.hh
               dginversemass.hh
               dgpass.hh
               insertfunction.hh
               insertfunctionpass.hh
               localdg.hh
               pass.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${pass_HEADERS} DESTINATION ${passdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
