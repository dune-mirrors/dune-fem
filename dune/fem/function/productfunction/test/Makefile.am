POLORD = 3
GRIDTYPE = YASPGRID
GRIDDIM = 2

LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(DUNEMPILIBS) $(LOCAL_LIBS)
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS) $(DUNEMPICPPFLAGS) -DPOLORDER=$(POLORD)

# remove this until this test is fixed again 
#ADAPTPROG = prod

check_PROGRAMS = $(ADAPTPROG)
# prod_SOURCES = productfunctiontest.cc

include $(top_srcdir)/am/global-rules
