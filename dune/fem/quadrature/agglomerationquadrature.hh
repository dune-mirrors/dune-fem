#ifndef DUNE_FEM_AGGLOMERATIONQUADRATURE_HH
#define DUNE_FEM_AGGLOMERATIONQUADRATURE_HH

#include <stack>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dune/fem/misc/threads/threadsafevalue.hh>

#include "quadrature.hh"
#include "elementquadrature.hh"

namespace Dune
{

  namespace Fem
  {


    /** \brief Agglomeration is a simple quadrature for polyhedral
     * cells based on sub triangulation  */
    template< typename GridPartImp, class IntegrationPointList >
    class Agglomeration
    {
    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the element quadrature
      enum { codimension = 0 };

      //! dimension of the grid
      enum { dimension = GridPartType :: dimension };

      //! dimension of the world
      enum { dimensionworld = GridPartType :: dimensionworld };

      //! type for reals (usually double)
      typedef typename GridPartType :: ctype RealType;

      typedef IntegrationPointList  IntegrationPointListType;

      typedef typename IntegrationPointListType :: CoordinateType  CoordinateType;

      typedef int QuadratureKeyType;

      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

    protected:
      typedef PolyhedronQuadrature< RealType, dimension > PolyhedronQuadratureType;
      typedef typename IntegrationPointListType :: IntegrationPointListType    IntegrationPointListImpl;

      typedef std::stack< std::unique_ptr< PolyhedronQuadratureType > >  PolyhedronQuadratureStorageType;

      // PolyhedronQuadrature storage
      static PolyhedronQuadratureStorageType& quadStorage()
      {
        static ThreadSafeValue< PolyhedronQuadratureStorageType > storage;
        return *storage;
      }

      // get object from stack or create new
      static PolyhedronQuadratureType* getObject( const GeometryType& type )
      {
        PolyhedronQuadratureStorageType& storage = quadStorage();
        if( storage.empty() )
        {
          return new PolyhedronQuadratureType( type, 0, IdProvider ::instance().newId() );
        }
        else
        {
          PolyhedronQuadratureType* quad = storage.top().release();
          assert( quad );
          storage.pop();
          return quad;
        }
      }

      // push object to stack or delete
      static void pushObject( const PolyhedronQuadratureType* quad )
      {
        PolyhedronQuadratureType* polyQuad = const_cast< PolyhedronQuadratureType* > (quad);
        PolyhedronQuadratureStorageType& storage = quadStorage();
        if( storage.size() < 20 )
        {
          storage.emplace( polyQuad );
        }
        else
        {
          delete polyQuad;
        }
      }

      // deleter object returning pointer to stack object
      struct Deleter
      {
        void operator ()(const IntegrationPointListImpl* quad)
        {
          const PolyhedronQuadratureType* polyQuad = static_cast< const PolyhedronQuadratureType* > ( quad );
          pushObject( polyQuad );
        }
      };

    public:
      //! returns quadrature points for polyhedral cells
      // this only works for 2d so far
      template <int dimw = dimensionworld >
      static std::enable_if_t<dimension == 2 && dimw == 2, IntegrationPointListType>
      computeQuadrature( const EntityType &entity, const QuadratureKeyType& quadKey )
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

        PolyhedronQuadratureType& quadImp = *(getObject( entity.type() ));

        quadImp.reset( order, subEntities * quadNop );

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

            quadImp.addQuadraturePoint( point, weight );
          }
        }
        // return a shared pointer with the correct deleter
        // removing the pointer to the stack
        std::shared_ptr< const IntegrationPointListImpl > quadPtr( &quadImp, Deleter() );
        return IntegrationPointListType( quadPtr );
      }

      template <int dimw = dimensionworld >
      static std::enable_if_t<dimension != 2 || dimw != 2, IntegrationPointListType>
      computeQuadrature( const EntityType &entity, const QuadratureKeyType& quadKey )
      {
        assert(false);
        return IntegrationPointListType( {} );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
