# don't follow the full GNU-standard; we definitely need automake >= 1.5
AUTOMAKE_OPTIONS = foreign 1.5

SUBDIRS = ode oemsolver test

noinst_LTLIBRARIES = libsolver.la

# construct only from sub-libraries
libsolver_la_SOURCES =
# if par DG was found, then ODE solver are not included  
if PARDG 
libsolver_la_LIBADD =
else 
libsolver_la_LIBADD = ode/libodesolver.la
endif

solverdir =  $(includedir)/dune/fem/solver
solver_HEADERS = inverseoperators.hh \
                 istlsolver.hh \
		 odesolver.hh \
		 oemsolver.hh \
		 multistep.hh \
		 pardg.hh \
                 rungekutta.hh \
		 timeprovider.hh
