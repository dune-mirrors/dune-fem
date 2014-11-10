dune_install(cachingpointlist.hh cachingquadrature.hh cornerpointset.hh
             defaultquadratures.hh dunequadratures.hh elementpointlist.hh
             elementpointlistbase.hh elementquadrature.hh femquadratures.hh
             femquadratures_inline.hh gausspoints.hh gausspoints_implementation.hh
             idprovider.hh integrator.hh quadrature.hh
             quadratureimp.hh quadratureimp_inline.hh pyramidpoints.hh
             quadprovider.hh simplexpoints.hh intersectionquadrature.hh)

dune_add_subdirs(caching test)

dune_add_library(quadrature OBJECT 
    caching/twistutility.cc 
    caching/topology.cc )
