# IGNORE_HEADERS

# linker flags
LDADD = $(ALL_PKG_LDFLAGS) $(ALL_PKG_LIBS) $(LOCAL_LIBS)
# preprocessor flags
AM_CPPFLAGS = $(ALL_PKG_CPPFLAGS)
AM_CXXFLAGS = -O0 -g -Wall $(ALL_PKG_CXXFLAGS)

TESTS = test-streams
check_PROGRAMS = test-streams

test_streams_SOURCES = test-streams.cc

streamstestdir =  $(includedir)/dune/fem/io/streams/test
streamstest_HEADERS = test-streams.cc 
