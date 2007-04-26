noinst_LTLIBRARIES = libodesolver.la
libodesolver_la_CXXFLAGS = @AM_CPPFLAGS@ 
libodesolver_la_LIBADD  =
libodesolver_la_SOURCES = odeobjectfiles.cc 

odedir = $(includedir)/dune/fem/solver/ode
ode_HEADERS = blas.hpp  iterative_solver.hpp limiter.hpp \
quasi_exact_newton.hpp linear_solver.hpp random.hpp \
communicator.hpp  matrix.hpp subblas.hpp \
dynamical_object.hpp   newton.hpp  thread.hpp \
emptycommunicator.hpp  nonlinear_solver.hpp timer.hpp \
function.hpp ode_solver.hpp  vector.hpp

include $(top_srcdir)/am/global-rules
