set(miscdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/misc)
set(SUBDIRS  petsc test threads)
set(misc_HEADERS  alugridwriter.hh
               bartonnackmaninterface.hh
               boundaryidentifier.hh
               boundaryidprovider.hh
               capabilities.hh
               checkgeomaffinity.hh
               codimmap.hh
               combineinterface.hh
               commhelper.hh
               debug.hh
               double.hh
               engineconcept.hh
               femeoc.hh
               femeoctable.hh
               femtimer.cc
               femtimer.hh
               field.hh
               fieldmatrixhelper.hh
               fmatrixconverter.hh
               functor.hh
               griddeclaration.hh
               gridname.hh
               gridobjectstreams.hh
               gridsolution.hh
               gridwidth.hh
               h1norm.hh
               l1error.hh
               l1norm.hh
               l2error.hh
               l2norm.hh
               linesegmentsampler.hh
               lpnorm.hh
               mapgeomtype.hh
               metaprogramming.hh
               mpimanager.hh
               nil.hh
               nonconformitylevel.hh
               referenceelementprovider.hh
               selectvalue.hh
               subobjects.hh
               suite.cc
               suite.hh
               test.cc
               test.hh
               threaditerator.hh
               threadmanager.hh
               validator.hh)

dune_add_library(misc SOURCES femtimer.cc test.cc suite.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${misc_HEADERS} DESTINATION ${miscdir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
