set(SUBDIRS  common localdg test)
set(passdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/pass)
set(pass_HEADERS  applylocaloperator.hh
               createpass.hh
               dginversemass.hh
               insertfunction.hh
               localdg.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${pass_HEADERS} DESTINATION ${passdir})
dune_add_subdirs(${SUBDIRS})
