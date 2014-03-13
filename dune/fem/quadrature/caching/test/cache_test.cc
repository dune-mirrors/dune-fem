#include <config.h>

#include <dune/grid/sgrid.hh>

#include <dune/fem/gridpart/leafgridpart.hh>

#ifdef ENABLE_ALUGRID
#include <dune/grid/alugrid.hh>
#endif

#ifdef ENABLE_ALBERTA
#include <dune/grid/albertagrid.hh>
#endif

#ifdef ENABLE_UGGRID
#define USE_PRISMGRID
#include <dune/grid/uggrid.hh>
#endif

#ifdef USE_PRISMGRID
//#include <dune/prismgrid/dgfgridtype.hh>
#endif

#include <dune/geometry/referenceelements.hh>

#include "cache_test.hh"

namespace Dune
{
  namespace Fem
  {

    void CacheProvider_Test::run ()
    {
      hexaTest();
      tetraTest();
      prismTest();
      triangleTest();
      quadTest();
    }

    void CacheProvider_Test::prismTest ()
    {
#ifdef USE_PRISMGRID
      const int dim = 3;
      const int codim = 1;

      //typedef PrismGrid<dim, dim> GridType;
      typedef UGGrid< dim > GridType;
      typedef Dune::Fem::LeafGridPart< GridType > GridPartType;

      typedef CacheProvider< GridPartType, codim > CacheProviderType;
      typedef Quadrature< GridPartType::ctype, GridPartType::dimension-1 > QuadratureType;

      typedef CacheProviderType::MapperType MapperType;
      typedef PointProvider< double, dim, codim > PointProviderType;
      typedef PointProviderType::GlobalPointVectorType PointVectorType;

      GeometryType elemGeo = GeometryType( GeometryType::prism, 3 );

      // Get reference element
      typedef Dune::ReferenceElement< double, dim > RefElement;
      const RefElement &refElement = Dune::ReferenceElements< double, dim >::general( elemGeo );

      // Loop over all faces
      for( int i = 0; i < refElement.size( codim ); ++i )
      {
        GeometryType faceGeo = (refElement.size( i, codim, dim ) == dim) ?
                               GeometryType( GeometryType::simplex, 2 ) :
                               GeometryType( GeometryType::cube, 2 );

        // Build quadrature
        QuadratureType quad( faceGeo, 3 );

        // Ask for one mapper so that the points get registered
        CacheProviderType::getMapper( quad, elemGeo, 0, 0 );

        const PointVectorType &points
          = PointProviderType::getPoints( quad.id(), elemGeo );

        const MapperType &m
          = CacheProviderType::getMapper( quad, elemGeo, i, 0 );

        _test( m.size() == (size_t) quad.nop());

        // Loop over all points
        const RefElement::Codim< codim >::Geometry refEmbedding = refElement.geometry< codim >( i );
        for( std::size_t j = 0; j < m.size(); ++j )
        {
          const FieldVector< double, dim > qpGlobal = refEmbedding.global( quad.point( j ) );
          for( int d = 0; d < dim; ++d )
            _floatTest( points[ m[ j ] ][ d ], qpGlobal[ d ] );
        }
      }
#endif
    }

    void CacheProvider_Test::hexaTest ()
    {
#ifdef ENABLE_ALUGRID
      const int dim = 3;
      const int codim = 1;

      typedef ALUCubeGrid< 3, 3 > GridType;
      typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
      typedef CacheProvider< GridPartType, codim > CacheProviderType;
      typedef Quadrature< GridPartType::ctype, GridPartType::dimension-1 > QuadratureType;
      //typedef CacheProviderType::QuadratureType QuadratureType;
      typedef CacheProviderType::MapperType MapperType;
      typedef PointProvider< double, dim, codim > PointProviderType;
      typedef PointProviderType::GlobalPointVectorType PointVectorType;

      GeometryType elemGeo = GeometryType( GeometryType::cube, 3 );
      GeometryType faceGeo = GeometryType( GeometryType::cube, 2 );

      // Get reference element
      typedef Dune::ReferenceElement< double, dim > RefElement;
      const RefElement &refElement = Dune::ReferenceElements< double, dim >::general( elemGeo );

      // Build quadrature
      QuadratureType quad( faceGeo, 3 );

      // Ask for one mapper so that the points get registered
      CacheProviderType::getMapper( quad, elemGeo, 0, 0 );

      const PointVectorType &points
        = PointProviderType::getPoints( quad.id(), elemGeo );

      // Loop over all faces
      for( int i = 0; i < refElement.size( codim ); ++i )
      {
        const MapperType &m
          = CacheProviderType::getMapper( quad, elemGeo, i, 0 );

        _test( m.size() == (size_t) quad.nop());
        // Loop over all points
        const RefElement::Codim< codim >::Geometry refEmbedding = refElement.geometry< codim >( i );
        for( std::size_t j = 0; j < m.size(); ++j )
        {
          const FieldVector< double, dim > qpGlobal = refEmbedding.global( quad.point( j ) );
          for( int d = 0; d < dim; ++d )
            _floatTest( points[ m[ j ] ][ d ], qpGlobal[ d ] );
        }
      }
#endif
    }

    void CacheProvider_Test::tetraTest ()
    {
#ifdef ENABLE_ALUGRID
      const int dim = 3;
      const int codim = 1;

      typedef ALUSimplexGrid< 3, 3 > GridType;
      typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
      typedef CacheProvider< GridPartType, codim > CacheProviderType;
      //typedef CacheProviderType::QuadratureType QuadratureType;
      typedef Quadrature< GridPartType::ctype, GridPartType::dimension-1 > QuadratureType;
      typedef CacheProviderType::MapperType MapperType;
      typedef PointProvider< double, dim, codim > PointProviderType;
      typedef PointProviderType::GlobalPointVectorType PointVectorType;

      GeometryType elemGeo = GeometryType( GeometryType::simplex, 3 );
      GeometryType faceGeo = GeometryType( GeometryType::simplex, 2 );

      // Get reference element
      typedef Dune::ReferenceElement< double, dim > RefElement;
      const RefElement &refElement = Dune::ReferenceElements< double, dim >::general( elemGeo );

      // Build quadrature
      QuadratureType quad( faceGeo, 3 );

      // Ask for one mapper so that the points get registered
      CacheProviderType::getMapper( quad, elemGeo, 0, 0 );

      const PointVectorType &points
        = PointProviderType::getPoints( quad.id(), elemGeo );

      // Loop over all faces
      for( int i = 0; i < refElement.size( codim ); ++i )
      {
        const MapperType &m = CacheProviderType::getMapper( quad, elemGeo, i, 0 );

        _test( m.size() == (size_t) quad.nop());
        // Loop over all points
        const RefElement::Codim< codim >::Geometry refEmbedding = refElement.geometry< codim >( i );
        for( std::size_t j = 0; j < m.size(); ++j )
        {
          const FieldVector< double, dim > qpGlobal = refEmbedding.global( quad.point( j ) );
          for( int d = 0; d < dim; ++d )
            _floatTest( points[ m[ j ] ][ d ], qpGlobal[ d ] );
        }
      }
#endif
    }

    void CacheProvider_Test::triangleTest ()
    {
#ifdef ENALBE_ALBERTA
      const int dim = 2;
      const int codim = 1;

      typedef AlbertaGrid< dim > GridType;
      typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
      typedef CacheProvider< GridPartType, codim > CacheProviderType;
      //typedef CacheProviderType::QuadratureType QuadratureType;
      typedef Quadrature< GridPartType::ctype, GridPartType::dimension-1 > QuadratureType;
      typedef CacheProviderType::MapperType MapperType;
      typedef PointProvider< double, dim, codim > PointProviderType;
      typedef PointProviderType::GlobalPointVectorType PointVectorType;

      GeometryType elemGeo = GeometryType( GeometryType::simplex, 2 );
      GeometryType faceGeo = GeometryType( GeometryType::simplex, 1 );

      // Get reference element
      typedef Dune::ReferenceElement< double, dim > RefElement;
      const RefElement &refElement = Dune::ReferenceElements< double, dim >::general( elemGeo );

      // Build quadrature
      QuadratureType quad( faceGeo, 3 );

      // Ask for one mapper so that the points get registered
      CacheProviderType::getMapper( quad, elemGeo, 0, 0 );

      const PointVectorType &points
        = PointProviderType::getPoints( quad.id(), elemGeo );

      // Loop over all faces
      for( int i = 0; i < refElement.size( codim ); ++i )
      {
        const MapperType &m = CacheProviderType::getMapper( quad, elemGeo, i, 0 );

        _test( m.size() == (size_t) quad.nop());
        // Loop over all points
        const RefElement::Codim< codim >::Geometry refEmbedding = refElement.geometry< codim >( i );
        for( size_t j = 0; j < m.size(); ++j )
        {
          const FieldVector< double, dim > qpGlobal = refEmbedding.global( quad.point( j ) );
          for( int d = 0; d < dim; ++d )
            _floatTest( points[ m[ j ] ][ d ], qpGlobal[ d ] );
        }
      }
#endif
    }

    void CacheProvider_Test::quadTest ()
    {
      const int dim = 2;
      const int codim = 1;

      typedef SGrid< dim, dim > GridType;
      typedef Dune::Fem::LeafGridPart< GridType > GridPartType;
      typedef CacheProvider< GridPartType, codim > CacheProviderType;
      typedef Quadrature< GridPartType::ctype, GridPartType::dimension-1 > QuadratureType;
      //typedef CacheProviderType::QuadratureType QuadratureType;
      typedef CacheProviderType::MapperType MapperType;
      typedef PointProvider< double, dim, codim > PointProviderType;
      typedef PointProviderType::GlobalPointVectorType PointVectorType;

      GeometryType elemGeo = GeometryType( GeometryType::cube, 2 );
      GeometryType faceGeo = GeometryType( GeometryType::cube, 1 );

      // Get reference element
      typedef Dune::ReferenceElement< double, dim > RefElement;
      const RefElement &refElement = Dune::ReferenceElements< double, dim >::general( elemGeo );

      // Build quadrature
      QuadratureType quad( faceGeo, 5 );

      // Ask for one mapper so that the points get registered
      CacheProviderType::getMapper( quad, elemGeo, 0, 0 );

      const PointVectorType &points
        = PointProviderType::getPoints( quad.id(), elemGeo );

      // Loop over all faces
      for( int i = 0; i < refElement.size( codim ); ++i )
      {
        const MapperType &m = CacheProviderType::getMapper( quad, elemGeo, i, 0 );

        _test( m.size() == (size_t) quad.nop());
        // Loop over all points
        const RefElement::Codim< codim >::Geometry refEmbedding = refElement.geometry< codim >( i );
        for( std::size_t j = 0; j < m.size(); ++j )
        {
          const FieldVector< double, dim > qpGlobal = refEmbedding.global( quad.point( j ) );
          for( int d = 0; d < dim; ++d )
            _floatTest( points[ m[ j ] ][ d ], qpGlobal[ d ] );
        }
      }
    }

  } // namespace Fem

} // namespace Dune
