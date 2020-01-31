#ifndef DUNE_FEM_AGGLOMERATIONQUADRATURE_HH
#define DUNE_FEM_AGGLOMERATIONQUADRATURE_HH

#include <dune/geometry/axisalignedcubegeometry.hh>

#include "quadrature.hh"
#include "elementquadrature.hh"

namespace Dune
{

  namespace Fem
  {


    /** \brief AgglomerationQuadrature is a simple quadrature for polyhedral
     * cells based on sub triangulation  */
    template< typename GridPartImp, class IntegrationPointList >
    class AgglomerationQuadrature
    {
    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the element quadrature
      enum { codimension = 0 };

      //! dimension of the world
      enum { dimension = GridPartType :: dimension };

      //! type for reals (usually double)
      typedef typename GridPartType :: ctype RealType;

      typedef IntegrationPointList  IntegrationPointListType;

      typedef typename IntegrationPointListType :: CoordinateType  CoordinateType;

      typedef int QuadratureKeyType;

      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

    public:
      static IntegrationPointListType computeQuadrature( const EntityType &entity, const QuadratureKeyType& quadKey )
      {
        typedef ElementQuadrature< GridPartImp, 0 > QuadratureType;
        Dune::GeometryType simplexType = Dune::GeometryTypes::simplex( dimension );

        typedef AxisAlignedCubeGeometry< RealType, dimension, dimension> CubeGeometryType;

        const auto& elemGeo = entity.geometry();

        // compute bounding box for setting up a cube quadrature
        CoordinateType lower = elemGeo.corner( 0 );
        CoordinateType upper = lower;
        const int corners = elemGeo.corners();
        for( int c = 1; c < corners; ++c )
        {
          const auto& corner = elemGeo.corner( c );
          for( int d=0; d< dimension; ++d )
          {
            lower[ d ] = std::min( lower[ d ], corner[ d ]);
            upper[ d ] = std::max( upper[ d ], corner[ d ]);
          }
        }

        //std::cout << "BBox = " << lower << " " << upper << std::endl;

        CubeGeometryType cubeGeom( lower, upper );

        QuadratureType quad( simplexType, quadKey );
        // TODO needs generalization
        const int subEntities = entity.subEntities( 1 );
        const int quadNop = quad.nop();
        int order = quad.order();

        std::vector< CoordinateType > points;
        std::vector< RealType > weights;

        points.reserve( subEntities * quadNop );
        weights.reserve( subEntities * quadNop );

        Dune::FieldMatrix<double,dimension,dimension> A( 0 ) ;
        const auto& center = entity.geometry().center();

        for( int i = 0; i<subEntities; ++i )
        {
          const auto subEntity = entity.template subEntity< 1 >( i );
          const auto& geom = subEntity.geometry();
          assert( geom.corners() == dimension );
          // setup transformation matrix, here setup transposed matrix
          for( int c = 0; c<dimension; ++c )
          {
            A[ c ]  = geom.corner( c );
            A[ c ] -= center;
          }

          // compute simplex volume / ref volume (which removes the 0.5)
          // abs is taken because determinant may be negative since we did not
          // care about the orientation
          double vol =  std::abs(A.determinant()) / elemGeo.volume();

          CoordinateType point;
          for( int qp = 0; qp < quadNop; ++qp )
          {
            // p = A^T * xHat + p_0 (A is stored as transposed)
            A.mtv( quad.point( qp ), point );
            point += center;

            point = cubeGeom.local( point );

            // scale weights with number of sub-triangles
            double weight = quad.weight( qp ) * vol;

            points.push_back( point );
            weights.push_back( weight );
          }
        }
        static PolyhedronQuadrature< RealType, dimension > quadImp(entity.type(), 0, IdProvider ::instance().newId() );
        quadImp.setQuadraturePoints( order, std::move(points), std::move( weights ) );
        return IntegrationPointListType( quadImp );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
