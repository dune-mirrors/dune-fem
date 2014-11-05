set(SUBDIRS  test)
set(filedir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/io/file)
set(file_HEADERS asciiparser.hh
                 dataoutput.hh
                 datawriter.hh
                 iointerface.hh
                 iolock.hh
                 iotuple.hh
                 latextablewriter.hh
                 persistencemanager.hh
                 vtkio.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${file_HEADERS} DESTINATION ${filedir})
dune_add_subdirs(${SUBDIRS})
