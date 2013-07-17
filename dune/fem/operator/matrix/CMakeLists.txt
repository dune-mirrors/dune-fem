set(matrixdir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator/matrix)
set(matrix_HEADERS  blockmatrix.hh columnobject.hh densematrix.hh  istlmatrix.hh
                 matrixsingleton.hh ontheflymatrix.hh
                 preconditionerwrapper.hh spmatrix.cc spmatrix.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${matrix_HEADERS} DESTINATION ${matrixdir})