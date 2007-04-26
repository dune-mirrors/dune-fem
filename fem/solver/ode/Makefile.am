noinst_LTLIBRARIES = libodesolver.la
libodesolver_la_CXXFLAGS = @AM_CPPFLAGS@ 
libodesolver_la_LIBADD  =
libodesolver_la_SOURCES = odeobjectfiles.cc 

CXXFLAGS=-O3 -DNDEBUG -funroll-loops -finline-functions --param \
max-inline-insns-single=4000 --param large-function-growth=4500 --param \
inline-unit-growth=4000 -ffast-math -fomit-frame-pointer -march=opteron -msse3 -mfpmath=sse

odedir = $(includedir)/dune/fem/solver/ode
ode_HEADERS = blas.hpp  iterative_solver.hpp limiter.hpp \
quasi_exact_newton.hpp linear_solver.hpp random.hpp \
communicator.hpp  matrix.hpp subblas.hpp \
dynamical_object.hpp   newton.hpp  thread.hpp \
emptycommunicator.hpp  nonlinear_solver.hpp timer.hpp \
function.hpp ode_solver.hpp  vector.hpp

include $(top_srcdir)/am/global-rules
