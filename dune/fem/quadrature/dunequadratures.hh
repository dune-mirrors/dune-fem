#ifndef DUNE_FEM_DUNEQUADRATURES_HH
#define DUNE_FEM_DUNEQUADRATURES_HH

//- Dune includes
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/fem/quadrature/quadratureimp.hh>

namespace Dune
{

  namespace Fem
  {

    /** \class QuadratureRulesFactory
     *  \ingroup Quadrature
     *  \brief quadrature implementation based on the standard DUNE quadratures
     *
     *  Though a factory by name, this is a quadrature implementation using the
     *  standard quadratures from DUNE grid to generate a list of quadrature
     *  points.
     */
    template< typename FieldImp, int dim >
    class QuadratureRulesFactory
    : public QuadratureImp< FieldImp, dim >
    {
    public:
      typedef FieldImp FieldType;

      enum { dimension = dim };

    private:
      typedef QuadratureRulesFactory< FieldType, dimension > ThisType;
      typedef QuadratureImp< FieldType, dimension > BaseType;

    protected:
      using BaseType :: addQuadraturePoint;

    public:
      typedef typename BaseType :: CoordinateType CoordinateType;

    protected:
      typedef QuadratureRule< FieldType, dimension > DuneQuadratureRuleType;

      //enum { highest_order_cube    = CubeQuadratureRule<ct,dim>::highest_order };
      //enum { highest_order_simplex = SimplexQuadratureRule<ct,dim>::highest_order };

      enum { highest_order = 44 };
      //(highest_order_cube < highest_order_simplex) ?
      //            highest_order_cube : highest_order_simplex };

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
      QuadratureRulesFactory( const GeometryType &geometry,
                              const int order,
                              const size_t id )
      : BaseType( id ),
        elementGeometry_( geometry )
      {
        // get gauss quadrature
        const DuneQuadratureRuleType &rule
          = QuadratureRules< FieldType, dimension >
            //:: rule( geometry, order, QuadratureType :: GaussLegendre );
            :: rule( geometry, order, QuadratureType :: GaussLobatto );

        order_ = rule.order();
        assert( order <= order_ );

        typedef typename DuneQuadratureRuleType :: iterator IteratorType;
        const IteratorType endit = rule.end();
        for( IteratorType it = rule.begin(); it != endit; ++it )
          addQuadraturePoint( (*it).position(), (*it).weight() );
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
        return highest_order;
      }
    };



    template< class FieldType, int dim >
    struct DuneQuadratureTraits
    {
      typedef QuadratureRulesFactory< FieldType, dim > SimplexQuadratureType;
      typedef QuadratureRulesFactory< FieldType, dim > CubeQuadratureType;

      typedef QuadratureImp< FieldType, dim > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };

    template< class FieldType >
    struct DuneQuadratureTraits< FieldType, 0 >
    {
      typedef QuadratureRulesFactory< FieldType, 0 > PointQuadratureType;

      typedef QuadratureImp< FieldType, 0 > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };

    template< class FieldType >
    struct DuneQuadratureTraits< FieldType, 1 >
    {
      typedef QuadratureRulesFactory< FieldType, 1 > LineQuadratureType;

      typedef QuadratureImp< FieldType, 1 > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };

    template< class FieldType >
    struct DuneQuadratureTraits< FieldType, 3 >
    {
      typedef QuadratureRulesFactory< FieldType, 3 > SimplexQuadratureType;
      typedef QuadratureRulesFactory< FieldType, 3 > CubeQuadratureType;

      typedef QuadratureRulesFactory< FieldType, 3 > PrismQuadratureType;
      typedef QuadratureRulesFactory< FieldType, 3 > PyramidQuadratureType;

      typedef QuadratureImp< FieldType, 3 > IntegrationPointListType;

      typedef int QuadratureKeyType ;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DUNEQUADRATURES_HH
