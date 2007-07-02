noinst_LTLIBRARIES = libodesolver.la

if MPI
libodesolver_la_CXXFLAGS = @AM_CPPFLAGS@ $(MPI_CPPFLAGS)
libodesolver_la_LIBADD  = $(MPI_LDFLAGS) $(MPILIBS)
else 
libodesolver_la_CXXFLAGS = @AM_CPPFLAGS@
libodesolver_la_LIBADD  = 
endif

libodesolver_la_SOURCES = odeobjectfiles.cc 

#CXXFLAGS= -O3 -DNDEBUG -funroll-loops -finline-functions --param \
#max-inline-insns-single=4000 --param large-function-growth=4500 --param \
#inline-unit-growth=4000 -ffast-math -march=opteron -msse3 -mfpmath=sse

odedir = $(includedir)/dune/fem/solver/ode
ode_HEADERS = blas.hpp  iterative_solver.hpp limiter.hpp \
quasi_exact_newton.hpp linear_solver.hpp random.hpp \
communicator.hpp mpicommunicator.hpp  matrix.hpp subblas.hpp \
dynamical_object.hpp   newton.hpp  thread.hpp \
emptycommunicator.hpp  nonlinear_solver.hpp timer.hpp \
function.hpp ode_solver.hpp  vector.hpp bicgstab.hpp 

EXTRA_DIST = bicgstab.cpp  dirk.cpp  matrix.cpp  quasi_exact_newton.cpp \
buffer.cpp          fgmres.cpp            newton.cpp      runge_kutta.cpp \
bulirsch_stoer.cpp  gmres.cpp             ode_solver.cpp  sirk.cpp \
cg.cpp              iterative_solver.cpp  pc_limiter.cpp  ssp.cpp \
communicator.cpp    lu_solver.cpp         qr_solver.cpp   vector.cpp 

include $(top_srcdir)/am/global-rules
