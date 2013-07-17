# do everything in the subdirs
set(SUBDIRS  common pass space function gridpart quadrature operator petsc solver storage io misc test)
set(noinst_LTLIBRARIES  libfem.la)
# construct only from sub-libraries
add_executable(libfem_la )
set(libfem_la_LIBADD  io/libio.la misc/libmisc.la space/libspace.la solver/libsolver.la quadrature/libquadrature.la)
# install all files to dune/fem
set(femdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem)
set(fem_HEADERS  version.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${fem_HEADERS} DESTINATION ${femdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
