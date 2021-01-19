# add custom target to build tool chain for python tests
dune_python_add_test(NAME dg_python
                     SCRIPT discontinuousgalerkin.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
