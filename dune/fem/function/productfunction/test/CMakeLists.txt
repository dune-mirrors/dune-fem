set(POLORD  3)
set(GRIDTYPE  YASPGRID)
set(GRIDDIM  2)
set(LDADD  ${ALL_PKG_LDFLAGS} ${ALL_PKG_LIBS} ${DUNEMPILIBS} ${LOCAL_LIBS})
set(AM_CPPFLAGS  ${ALL_PKG_CPPFLAGS} ${DUNEMPICPPFLAGS} -DPOLORDER=${POLORD})
# remove this until this test is fixed again
#ADAPTPROG = prod
set(check_PROGRAMS  ${ADAPTPROG})
add_executable(prod productfunctiontest.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
