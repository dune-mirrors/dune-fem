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
      template< typename FieldImp, int dim, class PointSet >
      class InterpolationQuadratureFactory
      : public Dune::Fem::QuadratureImp< FieldImp, dim >
      {
      public:
        typedef FieldImp FieldType;
        typedef PointSet PointSetType;

        enum { dimension = dim };

      private:
        typedef InterpolationQuadratureFactory< FieldType, dimension, PointSetType > ThisType;
        typedef Dune::Fem::QuadratureImp< FieldType, dimension > BaseType;

      protected:
        using BaseType :: addQuadraturePoint;

      public:
        typedef typename BaseType :: CoordinateType CoordinateType;

      protected:
        const GeometryType elementGeometry_;
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
          elementGeometry_( geometry )
        {
          // revert quadrature order to polynomial order
          const int polOrder = (order + 1)/2;
          PointSetType points( polOrder );
          if( geometry.isCube() )
          {
            points.buildCube();
          }
          else
          {
            DUNE_THROW(InvalidStateException,"InterpolationQuadratureFactory: Unsupported geometry type " << geometry);
          }

          order_ = points.order();

          for( unsigned int i=0; i<points.size(); ++i )
          {
            addQuadraturePoint( points[ i ].point(), points[ i ].weight() );
          }
        }

        /** \copydoc Dune::Fem::QuadratureImp::order
         */
        int order () const
        {
          return order_;
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
          return QuadratureRules<FieldType,dim>::
            maxOrder( Dune::GeometryTypes::cube(dim), Dune::QuadratureType::Enum(PointSet::pointSetId) );
        }
      };

      template< class FieldType, int dim, class PointSet >
      struct InterpolationQuadratureTraitsImpl
      {
        static const int pointSetId = PointSet::pointSetId;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > SimplexQuadratureType;
        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > CubeQuadratureType;

        typedef Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };

      template< class FieldType, class PointSet >
      struct InterpolationQuadratureTraitsImpl< FieldType, 0, PointSet >
      {
        static const int dim = 0;

        static const int pointSetId = PointSet::pointSetId;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > PointQuadratureType;

        typedef  Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };

      template< class FieldType, class PointSet >
      struct InterpolationQuadratureTraitsImpl< FieldType, 1, PointSet >
      {
        static const int dim = 1;

        static const int pointSetId = PointSet::pointSetId;

        typedef InterpolationQuadratureFactory< FieldType, dim, PointSet > LineQuadratureType;

        typedef  Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;

        typedef int QuadratureKeyType ;
      };

      template< class FieldType, class PointSet >
      struct InterpolationQuadratureTraitsImpl< FieldType, 3, PointSet >
      {
        static const int dim = 3;

        static const int pointSetId = PointSet::pointSetId;

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
  //  GaussLobatto point set (same quadrature as in dune-geometry
  //  but different point ordering)
  //
  //////////////////////////////////////////////////////////////////
  template <class  FieldType, int dim >
  using GaussLobattoQuadratureTraits = detail::InterpolationQuadratureTraitsImpl< FieldType, dim, GaussLobattoPointSet< FieldType, dim > >;

  //////////////////////////////////////////////////////////////////
  //
  //  GaussLegendre point set (same quadrature as in dune-geometry
  //  but different point ordering)
  //
  //////////////////////////////////////////////////////////////////
  template <class  FieldType, int dim >
  using GaussLegendreQuadratureTraits = detail::InterpolationQuadratureTraitsImpl< FieldType, dim, GaussLegendrePointSet< FieldType, dim > >;
#endif

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_LAGRANGEPOINTS_HH
