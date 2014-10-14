
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




message(AUTHOR_WARNING "TODO. Improve module test.")

