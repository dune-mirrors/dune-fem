# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILIBS) $(LOCAL_LIBS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS)
EXTRA_DIST = testparameterfile \
						 oldparameterfile \
						 param_DarcyStokes \
						 param_NavierStokes \
						 date.sh

TESTS = parametertest
check_PROGRAMS = parametertest

parametertest_SOURCES = parametertest.cc 


include $(top_srcdir)/am/global-rules

