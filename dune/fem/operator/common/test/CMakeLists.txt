set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} -DPOLORDER=${POLORDER})
set(check_PROGRAMS  mappings)
add_executable(mappings mappingtest.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
