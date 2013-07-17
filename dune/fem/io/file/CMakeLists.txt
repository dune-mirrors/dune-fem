set(SUBDIRS  test)
set(filedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/file)
set(file_HEADERS  asciiparser.hh
         datawriter.hh
         dataoutput.hh
         binarydataio.hh
         iointerface.hh
         iolock.hh
         iotuple.hh
         latextablewriter.hh
         persistencemanager.hh
         vtkio.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${file_HEADERS} DESTINATION ${filedir})
foreach(i ${SUBDIRS})
  if(${i} STREQUAL "test")
    set(opt EXCLUDE_FROM_ALL)
  endif(${i} STREQUAL "test")
  add_subdirectory(${i} ${opt})
endforeach(i ${SUBDIRS})
