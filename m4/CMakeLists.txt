set(ALLM4S compatibility.m4
           dune-fem.m4
           endian.m4
           papi.m4
           petsc.m4
           sionlib.m4
           thread_local_storage.m4
           xdr_uint64_t.m4)
set(aclocaldir  ${CMAKE_INSTALL_DATADIR}/aclocal)
set(aclocal_DATA  ${ALLM4S})
set(EXTRA_DIST  ${ALLM4S})
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${aclocal_DATA} DESTINATION ${aclocaldir})
