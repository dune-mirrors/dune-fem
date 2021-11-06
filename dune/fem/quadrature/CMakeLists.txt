dune_install(agglomerationquadrature.hh
             cachingpointlist.hh       cachingquadrature.hh cornerpointset.hh
             defaultquadratures.hh     dunequadratures.hh   elementpointlist.hh
             elementpointlistbase.hh   elementquadrature.hh femquadratures.hh
             femquadratures_inline.hh  gausspoints.hh       gausspoints_implementation.hh
             idprovider.hh             integrator.hh        interpolationquadrature.hh
             lumpingquadrature.hh      quadrature.hh        quadratureimp.hh
             quadratureimp_inline.hh   pyramidpoints.hh     quadprovider.hh
             intersectionquadrature.hh pardgsimplexquadrature.hh )

dune_add_subdirs(caching geometric test)

target_sources(dunefem PRIVATE
    pardgsimplexquadrature0d.cc
    pardgsimplexquadrature1d.cc
    pardgsimplexquadrature2d.cc
    pardgsimplexquadrature3d.cc
    caching/twistutility.cc
    caching/topology.cc
    )
