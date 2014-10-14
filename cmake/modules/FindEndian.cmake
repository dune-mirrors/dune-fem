# Module that checks whether Endian is available.
#
# Sets the follwing variable:
#
# HAVE_ENDIAN              True if Endian header available.
# SYSTEM_ENDIAN_HEADER     systems endian header in angle brackets.

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

#debug info (remove later)
#message( "Info: SYSTEM_ENDIAN_HEADER is set to: ${SYSTEM_ENDIAN_HEADER}")
