# searches for endian.h 
AC_DEFUN([DUNE_CHECK_ENDIANESS],[
  AC_CHECK_HEADERS([endian.h],[AC_DEFINE(SYSTEM_ENDIAN_HEADER,[<endian.h>],[ systems endian  header])])
  AC_CHECK_HEADERS([machine/endian.h],[AC_DEFINE(SYSTEM_ENDIAN_HEADER,[<machine/endian.h>],[ systems endian  header])])
])
