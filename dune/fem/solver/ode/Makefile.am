# install directory 
odedir = $(includedir)/dune/fem/solver/ode

# if parDG is used, don't include ode solver 
if PARDG
# only install these headers 
ode_HEADERS = communicator.hpp mpicommunicator.hpp emptycommunicator.hpp 
noinst_LTLIBRARIES = libodesolver.la
libodesolver_la_SOURCES =

else 
# otherwise install all 
noinst_LTLIBRARIES = libodesolver.la
if MPI
libodesolver_la_CXXFLAGS = @AM_CPPFLAGS@ $(MPI_CPPFLAGS)
libodesolver_la_LIBADD  = $(MPI_LDFLAGS) $(MPILIBS)
else 
libodesolver_la_CXXFLAGS = @AM_CPPFLAGS@
libodesolver_la_LIBADD  = 
endif

libodesolver_la_SOURCES = odeobjectfiles.cc 

ode_HEADERS = communicator.hpp mpicommunicator.hpp emptycommunicator.hpp \
blas.hpp  iterative_solver.hpp limiter.hpp \
quasi_exact_newton.hpp linear_solver.hpp random.hpp \
matrix.hpp subblas.hpp dynamical_object.hpp   newton.hpp  thread.hpp \
nonlinear_solver.hpp timer.hpp function.hpp ode_solver.hpp  vector.hpp \
quadrature.hpp 

EXTRA_DIST = communicator.cpp buffer.cpp \
bicgstab.cpp  dirk.cpp  matrix.cpp  quasi_exact_newton.cpp \
fgmres.cpp            newton.cpp      runge_kutta.cpp \
bulirsch_stoer.cpp  gmres.cpp             ode_solver.cpp  sirk.cpp \
cg.cpp              iterative_solver.cpp  pc_limiter.cpp  ssp.cpp \
lu_solver.cpp         qr_solver.cpp   vector.cpp \
quadrature0d.cpp  quadrature1d.cpp  quadrature2d.cpp  quadrature3d.cpp
endif
# end if PARDG 

include $(top_srcdir)/am/global-rules
