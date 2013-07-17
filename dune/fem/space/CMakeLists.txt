set(spacedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space)
set(space_HEADERS  combineddiscretefunctionspace.hh
                combinedspace.hh
                dgspace.hh
                discontinuousgalerkin.hh
                finitevolume.hh
                fourier.hh
                fvspace.hh
                lagrange.hh
                lagrangespace.hh
                padaptivespace.hh
                rannacherturek.hh
                reducedbasisspace.hh)
set(SUBDIRS  basefunctions basisfunctionset
          combinedspace common
          discontinuousgalerkin dofmapper
          finitevolume fourier lagrange mapper
          padaptivespace rannacherturek
          reducedbasisspace shapefunctionset
          test)
set(noinst_LTLIBRARIES  libspace.la)
# construct only from sub-libraries
add_executable(libspace_la )
set(libspace_la_LIBADD  discontinuousgalerkin/libdiscontinuousgalerkin.la
                     lagrange/liblagrange.la)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${space_HEADERS} DESTINATION ${spacedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
