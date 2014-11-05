# if parDG is used, don't include ode solver
if( PARDG)
    # only install these headers
    dune_install(communicator.hpp mpicommunicator.hpp emptycommunicator.hpp)
else( PARDG )
    # otherwise install all
    dune_install(communicator.hpp mpicommunicator.hpp emptycommunicator.hpp
                 blas.hpp iterative_solver.hpp limiter.hpp
                 quasi_exact_newton.hpp linear_solver.hpp random.hpp
                 matrix.hpp subblas.hpp dynamical_object.hpp
                 newton.hpp thread.hpp nonlinear_solver.hpp
                 timer.hpp function.hpp ode_solver.hpp 
                 vector.hpp quadrature.hpp)
endif( PARDG )
