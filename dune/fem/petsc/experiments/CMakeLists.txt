set(experimentsdir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/petsc/experiments)
set(SUBDIRS  petsc_matrix parallel_laplace commtest)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
