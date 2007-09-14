# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = ode oemsolver test

noinst_LTLIBRARIES = libsolver.la

# construct only from sub-libraries
libsolver_la_SOURCES =
libsolver_la_LIBADD = ode/libodesolver.la

solverdir =  $(includedir)/dune/fem/solver
solver_HEADERS =  istlsolver.hh rungekutta.hh odesolver.hh \
oemsolver.hh inverseoperators.hh pardg.hh 
