dune_install(blockdiagonal.hh istladapter.hh istloperator.hh
             petscoperator.hh spoperator.hh denserow.hh eigenoperator.hh
             numpyoperator.hh )

if( NOT dune-python_FOUND )
  exclude_from_headercheck( "numpyoperator.hh" )
endif()

dune_add_subdirs( test )
