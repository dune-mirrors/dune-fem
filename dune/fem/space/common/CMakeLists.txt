set(commondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/space/common)
set(common_HEADERS  adaptmanager.hh
     allgeomtypes.hh
     arrays.hh
     basesetlocalkeystorage.hh
     cachedcommmanager.hh
     capabilities.hh
     commindexmap.hh
     commoperations.hh
     communicationmanager.hh
     datacollector.hh
     defaultcommhandler.hh
     discretefunctionspace.hh
     dofmanager.hh
     dofstorage.hh
     functionspace.hh
     functionspaceinterface.hh
     interpolation.hh
     loadbalancer.hh
     localrestrictprolong.hh
     restrictprolongfunction.hh
     restrictprolonginterface.hh
     adaptcallbackhandle.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${common_HEADERS} DESTINATION ${commondir})