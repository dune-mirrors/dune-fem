# install directory
set(odedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/solver/ode)
# if parDG is used, don't include ode solver
if( PARDG)
    # only install these headers
    set(ode_HEADERS  communicator.hpp mpicommunicator.hpp emptycommunicator.hpp)
    
else( PARDG )
    # otherwise install all
    
    set(ode_HEADERS  communicator.hpp mpicommunicator.hpp emptycommunicator.hpp
    blas.hpp  iterative_solver.hpp limiter.hpp
    quasi_exact_newton.hpp linear_solver.hpp random.hpp
    matrix.hpp subblas.hpp dynamical_object.hpp   newton.hpp  thread.hpp
    nonlinear_solver.hpp timer.hpp function.hpp ode_solver.hpp  vector.hpp
    quadrature.hpp)
    set(EXTRA_DIST  communicator.cpp buffer.cpp
    bicgstab.cpp  dirk.cpp  matrix.cpp  quasi_exact_newton.cpp
    fgmres.cpp            newton.cpp      runge_kutta.cpp
    bulirsch_stoer.cpp  gmres.cpp             ode_solver.cpp  sirk.cpp
    cg.cpp              iterative_solver.cpp  pc_limiter.cpp  ssp.cpp
    lu_solver.cpp         qr_solver.cpp   vector.cpp
    quadrature0d.cpp  quadrature1d.cpp  quadrature2d.cpp  quadrature3d.cpp)
endif( PARDG )
# end if PARDG
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${ode_HEADERS} DESTINATION ${odedir})