# add custom target to build tool chain for python tests
dune_python_add_test(NAME advdiff_python
                     SCRIPT advectiondiffusion.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     TIMEOUT 2000)
