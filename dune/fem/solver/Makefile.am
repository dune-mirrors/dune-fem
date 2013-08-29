solverdir =  $(includedir)/dune/fem/solver
solver_HEADERS = cginverseoperator.hh \
                 diagonalpreconditioner.hh \
                 inverseoperators.hh \
                 istlsolver.hh \
                 odesolver.hh \
                 odesolverinterface.hh \
                 oemsolver.hh \
                 multistep.hh \
                 newtoninverseoperator.hh \
                 pardg.hh \
                 petscsolver.hh \
                 rungekutta.hh \
                 timeprovider.hh \
                 umfpacksolver.hh

SUBDIRS = ode oemsolver rungekutta test

noinst_LTLIBRARIES = libsolver.la

libsolver_la_SOURCES =
libsolver_la_LIBADD = rungekutta/librungekutta.la ode/libodesolver.la

include $(top_srcdir)/am/global-rules
