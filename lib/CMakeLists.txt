target_link_libraries(dunefem PUBLIC
  ${PAPI_LIBRARIES}
  ${SIONLIB_LIBRARIES}
  ${DUNE_LIBS})

add_dune_mpi_flags(dunefem)

install(FILES DESTINATION ${CMAKE_INSTALL_BINDIR}/../lib)
