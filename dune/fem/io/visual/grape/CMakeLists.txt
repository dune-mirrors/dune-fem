set(SUBDIRS  datadisp)
set(grapedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/visual/grape)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
dune_add_subdirs(${SUBDIRS})
