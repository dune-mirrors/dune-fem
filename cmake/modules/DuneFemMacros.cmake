
#make libdunefem known locally
set(LOCAL_LIBS "${PROJECT_BINARY_DIR}/lib/libdunefem.a"
		CACHE STRING "path to local libs in dune-fem" )

#find endian headers
check_include_file_cxx(endian.h HAVE_ENDIAN_HEADER_HH)
if( HAVE_ENDIAN_HEADER_HH )
  find_path(SYSTEM_ENDIAN_HEADER_PATH
    NAMES endian.h )
  set(SYSTEM_ENDIAN_HEADER "<endian.h>"
    CACHE STRING "Path where endian.h was found")
else(HAVE_ENDIAN_HEADER_HH)
  check_include_file_cxx(maschine/endian.h HAVE_ENDIAN_MACHINE_HEADER_HH)
  if( HAVE_ENDIAN_MACHINE_HEADER_HH )
    find_path(SYSTEM_ENDIAN_HEADER_PATH
      NAMES endian.h )
    set(SYSTEM_ENDIAN_HEADER "<machine/endian.h>"
      CACHE STRING "Path where endian.h was found")
  endif(HAVE_ENDIAN_MACHINE_HEADER_HH)
endif(HAVE_ENDIAN_HEADER_HH)
add_definitions("-DSYSTEM_ENDIAN_HEADER=${SYSTEM_ENDIAN_HEADER}")

#todo: there seems to be no cmake equivalent for
#xdr issues in dune-common. Thus, the following lines are only for completeness
#Please check whether it works and/or delete this stuff.
find_path(DUNE_PATH_XDR
  NAMES rpc/rpc.h )
if(DUNE_PATH_XDR)
  include(CheckFunctionExists)
  check_function_exists(xdr_uint64_t XDR_UINT64_FUNC_EXISTS)
  if( XDR_UINT64_FUNC_EXISTS )
    set(XDR_UINT64_FUNC "xdr_uint64_t"
      CACHE STRING "xdr_unit64_t routine")
  else(XDR_UINT64_FUNC_EXISTS )
    check_function_exists(xdr_u_int64_t XDR_U_INT64_FUNC_EXISTS)
    if( XDR_U_INT64_FUNC_EXISTS )
      set(XDR_UINT64_FUNC "xdr_u_int64_t"
        CACHE STRING "xdr_unit64_t routine")
    endif(XDR_U_INT64_FUNC_EXISTS )
  endif(XDR_UINT64_FUNC_EXISTS )
endif(DUNE_PATH_XDR)




message(AUTHOR_WARNING "TODO. Improve module test.")

