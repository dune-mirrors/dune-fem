#ifndef DUNE_FEMPY_FEMPYQUADRATURES_HH
#define DUNE_FEMPY_FEMPYQUADRATURES_HH

//- Dune includes
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>

#include <dune/fem/quadrature/quadratureimp.hh>

namespace Dune
{

  namespace FemPy
  {
    template <class Field, int dim>
    class QuadraturePointsRegistry
    {
      typedef QuadraturePointsRegistry< Field, dim > ThisType;

      typedef Dune::QuadratureRule< Field, dim > QuadratureRuleType;
      typedef Dune::QuadraturePoint< Field, dim > QuadraturePointType;
      typedef typename QuadraturePointType :: Vector  CoordinateType;

      std::vector< QuadratureRuleType > rules_;

      QuadraturePointsRegistry(int maxOrder = 100 ) : rules_( maxOrder ) {}

      QuadratureRuleType& getRule(const int order )
      {
        return rules_[ order ];
      }

      static ThisType& instance()
      {
        static ThisType reg;
        return reg;
      }

    public:
      DUNE_EXPORT static const QuadratureRuleType& quadratureRule( const int order )
      {
        return instance().rules_[ order ];
      }


      template <class Rules>
      DUNE_EXPORT
      static void registerQuadratureRule( const Rules& rules, const int order, const Dune::GeometryType& geometry )
      {
        QuadratureRuleType& rule = instance().getRule( order );

        pybind11::object pyrule = rules( geometry );
        pybind11::object pyPW = pyrule.attr("get")();
        auto pointsWeights = pyPW.template cast<
           std::pair<pybind11::array_t<double>,
                     pybind11::array_t<double>> >();

        // check shape here...
        for( std::size_t i = 0, sz = pointsWeights.second.size(); i < sz; ++i )
        {
          CoordinateType hatx(0);
          for (std::size_t c=0;c<CoordinateType::size();++c)
            hatx[c] = pointsWeights.first.at(c,i);
          double weight = pointsWeights.second.at( i );

          rule.push_back( QuadraturePointType( hatx, weight) );
        }
      }
    };


    /** \class FempyQuadratureRulesFactory
     *  \ingroup Quadrature
     *  \brief quadrature implementation based on the standard DUNE quadratures
     *
     *  Though a factory by name, this is a quadrature implementation using the
     *  standard quadratures from DUNE grid to generate a list of quadrature
     *  points.
     */
    template< typename FieldImp, int dim >
    class FempyQuadratureRulesFactory
    : public Dune::Fem::QuadratureImp< FieldImp, dim >
    {
    public:
      typedef FieldImp FieldType;

      enum { dimension = dim };

    private:
      typedef FempyQuadratureRulesFactory< FieldType, dimension > ThisType;
      typedef Dune::Fem::QuadratureImp< FieldType, dimension > BaseType;

    protected:
      using BaseType :: addQuadraturePoint;

    public:
      typedef typename BaseType :: CoordinateType CoordinateType;

    protected:
      typedef Dune::QuadratureRule< FieldType, dimension > DuneQuadratureRuleType;
      typedef QuadraturePointsRegistry< FieldType, dimension >
        QuadratureRuleRegistryType;

      //enum { highest_order_cube    = CubeQuadratureRule<ct,dim>::highest_order };
      //enum { highest_order_simplex = SimplexQuadratureRule<ct,dim>::highest_order };

      enum { highest_order = 100 };
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
      FempyQuadratureRulesFactory( const GeometryType &geometry,
                                   const int order,
                                   const size_t id )
      : BaseType( id ),
        elementGeometry_( geometry )
      {
        const DuneQuadratureRuleType &rule =
          QuadratureRuleRegistryType::quadratureRule( order );

        order_ = order; // ???? rule.order();
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
    struct FempyQuadratureTraits
    {
      typedef FempyQuadratureRulesFactory< FieldType, dim > SimplexQuadratureType;
      typedef FempyQuadratureRulesFactory< FieldType, dim > CubeQuadratureType;

      typedef Dune::Fem::QuadratureImp< FieldType, dim > IntegrationPointListType;
    };

    template< class FieldType >
    struct FempyQuadratureTraits< FieldType, 0 >
    {
      typedef FempyQuadratureRulesFactory< FieldType, 0 > PointQuadratureType;

      typedef Dune::Fem::QuadratureImp< FieldType, 0 > IntegrationPointListType;
    };

    template< class FieldType >
    struct FempyQuadratureTraits< FieldType, 1 >
    {
      typedef FempyQuadratureRulesFactory< FieldType, 1 > LineQuadratureType;

      typedef Dune::Fem::QuadratureImp< FieldType, 1 > IntegrationPointListType;
    };

    template< class FieldType >
    struct FempyQuadratureTraits< FieldType, 3 >
    {
      typedef FempyQuadratureRulesFactory< FieldType, 3 > SimplexQuadratureType;
      typedef FempyQuadratureRulesFactory< FieldType, 3 > CubeQuadratureType;

      typedef FempyQuadratureRulesFactory< FieldType, 3 > PrismQuadratureType;
      typedef FempyQuadratureRulesFactory< FieldType, 3 > PyramidQuadratureType;

      typedef Dune::Fem::QuadratureImp< FieldType, 3 > IntegrationPointListType;
    };

  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FEMPYQUADRATURES_HH
