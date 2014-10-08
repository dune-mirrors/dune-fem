# install directory
set(odedir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/solver/ode)
# if parDG is used, don't include ode solver
if( PARDG)
    # only install these headers
    set(ode_HEADERS  communicator.hpp 
                     mpicommunicator.hpp 
                     emptycommunicator.hpp)
    
else( PARDG )
    # otherwise install all
    
    set(ode_HEADERS communicator.hpp 
                    mpicommunicator.hpp 
                    emptycommunicator.hpp
                    blas.hpp 
                    iterative_solver.hpp
                    limiter.hpp
                    quasi_exact_newton.hpp
                    linear_solver.hpp 
                    random.hpp
                    matrix.hpp
                    subblas.hpp
                    dynamical_object.hpp
                    newton.hpp 
                    thread.hpp
                    nonlinear_solver.hpp
                    timer.hpp 
                    function.hpp
                    ode_solver.hpp 
                    vector.hpp
                    quadrature.hpp)
endif( PARDG )
# end if PARDG
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${ode_HEADERS} DESTINATION ${odedir})
