#ifndef DUNE_FEM_INTERPOLATIONQUADRATURE_HH
#define DUNE_FEM_INTERPOLATIONQUADRATURE_HH

//- Dune includes
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/fem/quadrature/quadratureimp.hh>
#include <dune/fem/space/localfiniteelement/quadratureinterpolation.hh>

namespace Dune
{

  namespace Fem
  {

    namespace detail
    {
      /** \class InterpolationQuadratureFactory
       *  \ingroup Quadrature
       *  \brief
       *
       */
      template< typename FieldImp, int dim, template <class,unsigned int> class PointSet >
      class InterpolationQuadratureFactory
      : public Dune::Fem::QuadratureImp< FieldImp, dim >
      {
      public:
        typedef FieldImp FieldType;
        typedef PointSet< FieldType, dim > PointSetType;

        enum { dimension = dim };

      private:
        typedef InterpolationQuadratureFactory< FieldType, dimension, PointSet > ThisType;
        typedef Dune::Fem::QuadratureImp< FieldType, dimension > BaseType;

      protected:
        typedef typename BaseType :: ElementCoordinateType ElementCoordinateType;
        using BaseType :: addQuadraturePoint;

      public:
        typedef typename BaseType :: CoordinateType CoordinateType;

      protected:
        const GeometryType elementGeometry_;
        mutable size_t numElementInterpolPoints_;
        int order_;

      public:
        /** \brief constructor filling the list of points and weights
         *
         *  \param[in]  geometry  geometry type for which a quadrature is desired
         *  \param[in]  order     desired order (provided by the user)
         *  \param[in]  id        unique identifier (provided by QuadratureProvider)
         */
        InterpolationQuadratureFactory( const GeometryType &geometry,
                                        const int order,
                                        const size_t id )
        : BaseType( id ),
          elementGeometry_( geometry ),
          numElementInterpolPoints_( 0 ) // this is only set when interpolationPoints are requested
        {
          // revert quadrature order to polynomial order
          if( geometry.isCube() )
          {
            auto points = PointSetType::buildCubeQuadrature( order );
            order_ = points.quadOrder();

            assert(order_ >= order);

            for( unsigned int i=0; i<points.size(); ++i )
            {
              addQuadraturePoint( points[ i ].point(), points[ i ].weight() );
            }
          }
          else
          {
            DUNE_THROW(InvalidStateException,"InterpolationQuadratureFactory: Unsupported geometry type " << geometry);
          }
        }

        /** \copydoc Dune::Fem::QuadratureImp::order
         */
        int order () const
        {
          return order_;
        }

        /** \copydoc Dune::Fem::QuadratureImp::interpolationPoints
         */
        virtual std::vector< ElementCoordinateType > interpolationPoints(const int reqDim ) const
        {
          // if requested dimension matches potential element dimension
          if( reqDim == (dim+1) && dim < 3 )
          {
            typedef PointSet< FieldType, dim+1 > ElementPointSet;
            auto points = ElementPointSet::buildCubeQuadrature( order_ );
            numElementInterpolPoints_ = points.size();
            std::vector< ElementCoordinateType > pts( numElementInterpolPoints_ );
            for( size_t i=0; i<numElementInterpolPoints_; ++i )
              pts[ i ] = points[ i ].point();
            return pts;
          }
          else
            return std::vector< ElementCoordinateType >();
        }

        /** \copydoc Dune::Fem::QuadratureImp::isFaceInterpolationQuadrature */
        virtual bool isFaceInterpolationQuadrature( const size_t numShapeFunctions ) const
        {
          // when numElementInterpolPoints_ is set this means we have a face
          // quadrature and then this is an interpolation quadrature if the
          // number of shape functions matches the number of element
          // interpolation points
          return numShapeFunctions == numElementInterpolPoints_;
        }

        /** \copydoc Dune::Fem::QuadratureImp::geometry
         */
        GeometryType geometryType () const
        {
          return elementGeometry_;
        }

        /** \brief maximal order of available quadratures
         */
        static unsigned int maxOrder ()
        {
          // if the quadrature is from dune-geometry use the info from there
          if constexpr( PointSetType::pointSetId < Dune::QuadratureType::size )
          {
            return QuadratureRules<FieldType,dim>::
              maxOrder( Dune::GeometryTypes::cube(dim), Dune::QuadratureType::Enum(PointSetType::pointSetId) );
          }
          else
          {
            return 20;
          }
        }
      };

      template< class FieldType, int dim, template <class,unsigned int> class PointSet >
      struct InterpolationQuadratureTraitsImpl
      {
        static const int pointSetId = PointSet<FieldType,dim>::pointSetId;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > SimplexQuadratureType;
        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > CubeQuadratureType;

        typedef Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };

      template< class FieldType, template <class,unsigned int> class PointSet >
      struct InterpolationQuadratureTraitsImpl< FieldType, 0, PointSet >
      {
        static const int dim = 0;

        static const int pointSetId = PointSet<FieldType,dim>::pointSetId;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > PointQuadratureType;

        typedef  Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };

      template< class FieldType, template <class,unsigned int> class PointSet >
      struct InterpolationQuadratureTraitsImpl< FieldType, 1, PointSet >
      {
        static const int dim = 1;

        static const int pointSetId = PointSet<FieldType,dim>::pointSetId;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > LineQuadratureType;

        typedef  Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };

      template< class FieldType, template <class,unsigned int> class PointSet >
      struct InterpolationQuadratureTraitsImpl< FieldType, 3, PointSet >
      {
        static const int dim = 3;

        static const int pointSetId = PointSet<FieldType,dim>::pointSetId;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > SimplexQuadratureType;
        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > CubeQuadratureType;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > PrismQuadratureType;
        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > PyramidQuadratureType;

        typedef  Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };
  } // end namespace detail

#if HAVE_DUNE_LOCALFUNCTIONS
  //////////////////////////////////////////////////////////////////
  //
  //  Equidistant point set
  //
  //////////////////////////////////////////////////////////////////
  template <class  FieldType, int dim >
  using EquidistantQuadratureTraits = detail::InterpolationQuadratureTraitsImpl< FieldType, dim, EquidistantPointSetDerived >;

  //////////////////////////////////////////////////////////////////
  //
  //  GaussLobatto point set (same quadrature as in dune-geometry
  //  but different point ordering)
  //
  //////////////////////////////////////////////////////////////////
  template <class  FieldType, int dim >
  using GaussLobattoQuadratureTraits = detail::InterpolationQuadratureTraitsImpl< FieldType, dim, GaussLobattoPointSet >;

  //////////////////////////////////////////////////////////////////
  //
  //  GaussLegendre point set (same quadrature as in dune-geometry
  //  but different point ordering)
  //
  //////////////////////////////////////////////////////////////////
  template <class  FieldType, int dim >
  using GaussLegendreQuadratureTraits = detail::InterpolationQuadratureTraitsImpl< FieldType, dim, GaussLegendrePointSet >;

  //////////////////////////////////////////////////////////////////
  //
  //  Cell centers point set
  //
  //////////////////////////////////////////////////////////////////
  template <class  FieldType, int dim >
  using CellCentersQuadratureTraits = detail::InterpolationQuadratureTraitsImpl< FieldType, dim, CellCentersPointSet >;

#endif

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGEPOINTS_HH
