set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS})
set(check_PROGRAMS  odetest odesystemtest)
add_executable(odetest odetest.cc)
add_executable(odesystemtest odesystemtest.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
