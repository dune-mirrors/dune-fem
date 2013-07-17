# IGNORE_HEADERS
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(check_PROGRAMS  tester)
add_executable(tester main.cc ../../misc/test.cc ../../misc/suite.cc
        storage_test.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
