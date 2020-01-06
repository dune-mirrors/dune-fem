#ifndef DUNE_FEM_AGGLOMERATIONQUADRATURE_HH
#define DUNE_FEM_AGGLOMERATIONQUADRATURE_HH

#include <dune/geometry/axisalignedcubegeometry.hh>

#include "quadrature.hh"
#include "elementquadrature.hh"

namespace Dune
{

  namespace Fem
  {

    /*! \class \docme
     */
    template< typename GridPartImp, int codim, template <class, int> class QuadratureTraits = DefaultQuadratureTraits >
    class AgglomerationQuadrature;


    /** \copydoc AgglomerationQuadrature */
    template< typename GridPartImp, template <class, int> class QuadratureTraits >
    class AgglomerationQuadrature< GridPartImp, 0, QuadratureTraits >
      : public ElementQuadrature< GridPartImp, 0, QuadratureTraits >
    {
      typedef AgglomerationQuadrature< GridPartImp, 0 > ThisType;
      typedef ElementQuadrature< GridPartImp, 0, QuadratureTraits > BaseType;
      typedef typename BaseType :: IntegrationTraits   IntegrationTraits;

    public:
      //! type of the grid partition
      typedef GridPartImp GridPartType;

      //! codimension of the element quadrature
      enum { codimension = 0 };

      //! dimension of the world
      enum { dimension = GridPartType :: dimension };

      //! type for reals (usually double)
      typedef typename GridPartType :: ctype RealType;

      //! type for coordinates in the codim-0 reference element
      typedef typename IntegrationTraits :: CoordinateType CoordinateType;

      typedef typename BaseType :: IntegrationPointListType  IntegrationPointListType;

      //! type of quadrature identifier on user side (default is the order of quadrature)
      typedef int QuadratureKeyType;

      //! type of the quadrature point
      typedef QuadraturePointWrapper< ThisType > QuadraturePointWrapperType;
      //! type of iterator
      typedef QuadraturePointIterator< ThisType > IteratorType;

      // for compatibility
      typedef typename GridPartType::template Codim< 0 >::EntityType EntityType;

      template <class Geom>
      void computeBoundingBox( const Geom& geom, CoordinateType& lower, CoordinateType& upper)
      {
        lower = geom.corner( 0 );
        upper = geom.corner( 0 );
        const int corners = geom.corners();
        for( int c = 1; c < corners; ++c )
        {
          const auto& corner = geom.corner( c );
          for( int d=0; d< dimension; ++d )
          {
            lower[ d ] = std::min( lower[ d ], corner[ d ]);
            upper[ d ] = std::max( upper[ d ], corner[ d ]);
          }
        }
      }

    protected:
      using BaseType :: quad_;

    public:
      IntegrationPointListType createQuad( const EntityType &entity, const QuadratureKeyType& quadKey, const bool checkGeomType )
      {
        const GeometryType geomType = entity.type();
        if( checkGeomType && ! geomType.isNone() )
        {
          return IntegrationPointListType( geomType, quadKey );
        }
        else // compute weights and points based on sub-triangulation
        {
          typedef ElementQuadrature< GridPartImp, 0 > QuadratureType;
          Dune::GeometryType simplexType = Dune::GeometryTypes::simplex( dimension );

          typedef AxisAlignedCubeGeometry< RealType, dimension, dimension> CubeGeometryType;

          CoordinateType lower;
          CoordinateType upper;
          const auto& elemGeo = entity.geometry();
          computeBoundingBox( elemGeo, lower, upper );

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
          static PolyhedronQuadrature< RealType, dimension > quadImp(geomType, 0, IdProvider ::instance().newId() );
          quadImp.setQuadraturePoints( order, std::move(points), std::move( weights ) );
          return IntegrationPointListType( quadImp );
        }
      }

      /*! \brief constructor
       *
       *  \param[in]  entity    entity, on whose reference element the quadrature
       *                        lives
       *  \param[in]  quadKey   desired minimal order of the quadrature or other means of quadrature identification
       *  \param[in]  checkType if true entity's geometry type is checked
       *                        and if not none then a standard quadrature is used (default is true)
       */
      AgglomerationQuadrature( const EntityType &entity, const QuadratureKeyType& quadKey, const bool checkType = true )
        : BaseType( createQuad( entity, quadKey, checkType ) )
      {
      }

      /*! \brief constructor
       *
       *  \param[in]  type    geometry type, on whose reference element the quadrature
       *                      lives
       *  \param[in]  quadKey desired minimal order of the quadrature or other means of quadrature identification
       */
      AgglomerationQuadrature( const GeometryType &type, const QuadratureKeyType& quadKey )
      {
        DUNE_THROW(InvalidStateException,"AgglomerationQuadrature cannot be created with a geometry type. Needs list of elements");
      }

      QuadraturePointWrapperType operator[] ( std::size_t i ) const
      {
        return QuadraturePointWrapperType( *this, i );
      }

      IteratorType begin () const noexcept { return IteratorType( *this, 0 ); }
      IteratorType end () const noexcept { return IteratorType( *this, this->nop() ); }
    };

    template<class GridPart, class Entity>
    static inline auto agglomerationQuadrature(const GridPart& gridPart, const Entity& entity, unsigned quadOrder)
    {
      using Quadrature = Dune::Fem::AgglomerationQuadrature<GridPart, 0>;
      return Quadrature(entity, quadOrder);
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_ELEMENTQUADRATURE_HH
