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

# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${lagrange_HEADERS} DESTINATION ${lagrangedir})
dune_add_subdirs(${SUBDIRS})
