dune_python_add_test(NAME pybackrest1
                     SCRIPT backrest1.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)
dune_python_add_test(NAME pybackrest2
                     SCRIPT backrest2.py
                     WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                     LABELS quick)

if( ${FEM_TORTURE_TESTS} )
  dune_python_add_test(NAME filteredgridview
                       SCRIPT filteredgridview.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       LABELS quick)
  dune_python_add_test(NAME levelspaces
                       SCRIPT levelspaces.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       LABELS quick)
  dune_python_add_test(NAME pyassembly
                       SCRIPT test_assembly.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       LABELS quick)
  dune_python_add_test(NAME pyscheme
                       SCRIPT scheme_assemble.py
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                       TIMEOUT 2000
                       LABELS quick)
endif()

# which to add to test?
# fenics-heat.py       test-quad.py
# fenics-poisson.py    test-tupledf.py
# codegen.py           interpolate.py       testd2.py           testufl2vtk.py
# fempy-bcrsmatrix.py  localcontexttest.py  testnVector.py      testvecbackends.py
# fenics-elasticity.py mixed.py             testoperator.py
