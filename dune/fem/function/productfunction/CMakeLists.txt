set(SUBDIRS  test)
set(productfunctiondir   ${CMAKE_INSTALL_INCLUDEDIR}/dune/fem/function/productfunction)
set(productfunction_HEADERS  productfunction.hh  
                             productfunctionimp.cc)
# include not needed for CMake
# include $(top_srcdir)/am/global-rules
install(FILES ${productfunction_HEADERS} DESTINATION ${productfunctiondir})
dune_add_subdirs(${SUBDIRS})
