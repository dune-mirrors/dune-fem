set(modeldir  ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/operator/model)
set(model_HEADERS  linearellipticmodel.hh 
                   boundarymodel.hh 
                   diffusionmodel.hh 
                   discretemodels.hh)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${model_HEADERS} DESTINATION ${modeldir})
