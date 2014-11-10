#xdr
set(XDR_ROOT "" CACHE STRING "root of header xdr.h")
find_path(DUNE_PATH_XDR
  NAMES xdr.h
  PATHS ${XDR_ROOT} /usr /usr/include
  PATH_SUFFIXES "rpc"
  DOC "path to header xdr.h"
)

if(DUNE_PATH_XDR)
  include(CheckFunctionExists)
  check_function_exists(xdr_uint64_t XDR_UINT64_FUNC_EXISTS)
  if( XDR_UINT64_FUNC_EXISTS )
    set( XDR_UINT64_FUNC "xdr_uint64_t")
  else(XDR_UINT64_FUNC_EXISTS )
    check_function_exists(xdr_u_int64_t XDR_U_INT64_FUNC_EXISTS)
    if( XDR_U_INT64_FUNC_EXISTS )
      set( XDR_UINT64_FUNC "xdr_u_int64_t")
    endif(XDR_U_INT64_FUNC_EXISTS )
  endif(XDR_UINT64_FUNC_EXISTS )
else()
  message("can not find xdr.h")
endif(DUNE_PATH_XDR)
mark_as_advanced(XDR_ROOT DUNE_PATH_XDR)

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "XDR"
  DEFAULT_MSG
  DUNE_PATH_XDR
)
