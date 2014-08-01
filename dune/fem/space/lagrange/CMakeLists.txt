set(AUTOMAKE_OPTIONS  foreign 1.5)
set(SUBDIRS  test)
set(lagrangedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/lagrange)
set(lagrange_HEADERS  adaptmanager.hh
                   dofmappercode.hh
                   genericlagrangepoints.hh
                   shapefunctionset.hh
                   capabilities.hh
                   genericbasefunctions.hh
                   lagrangepoints.hh
                   space.hh
                   declaration.hh
                   genericgeometry.hh
                   restrictprolong.hh
                   storage.hh)
set(noinst_LTLIBRARIES  liblagrange.la)
dune_add_library(lagrange SOURCES shapefunctionset.cc)
target_link_libraries(lagrange ${DUNE_LIBS})

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${lagrange_HEADERS} DESTINATION ${lagrangedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
