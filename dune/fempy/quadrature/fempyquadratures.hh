#ifndef DUNE_FEMPY_FEMPYQUADRATURES_HH
#define DUNE_FEMPY_FEMPYQUADRATURES_HH

//- Dune includes
#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#if USING_DUNE_PYTHON
#include <dune/python/pybind11/pybind11.h>
#include <dune/python/pybind11/numpy.h>
#endif

#include <dune/fem/quadrature/quadprovider.hh>
#include <dune/fem/quadrature/quadratureimp.hh>
#include <dune/fem/quadrature/defaultquadratures.hh>
#include <dune/fem/storage/singleton.hh>

namespace Dune
{

  namespace FemPy
  {

    // forward declaration
    template< class DefaultQuadrature >
    class FempyQuadratureRulesFactory;

    template< class FieldType, int dim >
    struct FempyQuadratureTraits
    {
      typedef Dune::Fem::DefaultQuadratureTraits< FieldType, dim >
        DefaultQuadTraits;

      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits:: SimplexQuadratureType > SimplexQuadratureType;
      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits:: CubeQuadratureType    > CubeQuadratureType;

      typedef typename DefaultQuadTraits::IntegrationPointListType  IntegrationPointListType;
      typedef typename CubeQuadratureType :: QuadratureKeyType      QuadratureKeyType;
    };

    template< class FieldType >
    struct FempyQuadratureTraits< FieldType, 0 >
    {
      typedef Dune::Fem::DefaultQuadratureTraits< FieldType, 0 >
        DefaultQuadTraits;

      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits::PointQuadratureType > PointQuadratureType;
      typedef typename DefaultQuadTraits::IntegrationPointListType  IntegrationPointListType;
      typedef typename PointQuadratureType :: QuadratureKeyType     QuadratureKeyType;
    };

    template< class FieldType >
    struct FempyQuadratureTraits< FieldType, 1 >
    {
      typedef Dune::Fem::DefaultQuadratureTraits< FieldType, 1 >
        DefaultQuadTraits;

      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits::LineQuadratureType > LineQuadratureType;

      typedef typename DefaultQuadTraits::IntegrationPointListType IntegrationPointListType;
      typedef typename LineQuadratureType :: QuadratureKeyType     QuadratureKeyType;
    };

    template< class FieldType >
    struct FempyQuadratureTraits< FieldType, 3 >
    {
      typedef Dune::Fem::DefaultQuadratureTraits< FieldType, 3 >
        DefaultQuadTraits;

      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits:: SimplexQuadratureType > SimplexQuadratureType;
      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits:: CubeQuadratureType    > CubeQuadratureType;
      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits:: PrismQuadratureType   > PrismQuadratureType;
      typedef FempyQuadratureRulesFactory< typename DefaultQuadTraits:: PyramidQuadratureType > PyramidQuadratureType;

      typedef typename DefaultQuadTraits::IntegrationPointListType  IntegrationPointListType;
      typedef typename CubeQuadratureType :: QuadratureKeyType      QuadratureKeyType;
    };


    template <class Field, int dim>
    struct QuadraturePointsRegistry
    {
      typedef QuadraturePointsRegistry< Field, dim > ThisType;

      typedef Dune::QuadratureRule< Field, dim > QuadratureRuleType;
      typedef Dune::QuadraturePoint< Field, dim > QuadraturePointType;
      typedef typename QuadraturePointType :: Vector  CoordinateType;

    public:
      // defined in fem/quadrature/quadprovider.hh
      typedef Dune::Fem::FemQuadratureKey  QuadratureKeyType ;

    protected:
      std::map< QuadratureKeyType, QuadratureRuleType > rules_;

      QuadratureRuleType& getRule(const QuadratureKeyType& key)
      {
        return rules_[ key ];
      }

      bool quadExists( const QuadratureKeyType& key ) const
      {
        return rules_.find( key ) != rules_.end();
      }

      static ThisType& instance()
      {
        return Dune::Fem::Singleton< ThisType >::instance();
      }

    public:
      static const QuadratureRuleType& quadratureRule( const QuadratureKeyType& key )
      {
        assert( quadratureExists( key ) );
        return instance().rules_[ key ];
      }

      static bool quadratureExists( const QuadratureKeyType& key )
      {
        return instance().quadExists( key );
      }

      template <class Rules>
      static void registerQuadratureRule( const Rules& rules, const QuadratureKeyType& key, const Dune::GeometryType& geometry )
      {
#if USING_DUNE_PYTHON
        QuadratureRuleType& rule = instance().getRule( key );

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
#else
        DUNE_THROW(Dune::InvalidStateException,"registerQuadratureRule only works when USING_DUNE_PYTHON == 1");
#endif
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
    template< class DefaultQuadrature >
    class FempyQuadratureRulesFactory
    : public Dune::Fem::QuadratureImp< typename DefaultQuadrature::FieldType, DefaultQuadrature::dimension >
    {
    private:
      typedef FempyQuadratureRulesFactory< DefaultQuadrature > ThisType;
      typedef Dune::Fem::QuadratureImp< typename DefaultQuadrature::FieldType, DefaultQuadrature::dimension > BaseType;

    public:
      typedef typename  BaseType :: FieldType   FieldType;
      using BaseType :: dimension;

    protected:
      using BaseType :: addQuadraturePoint;

      typedef Dune::QuadratureRule< FieldType, dimension >      DuneQuadratureRuleType;
      typedef QuadraturePointsRegistry< FieldType, dimension >  QuadratureRuleRegistryType;

    public:
      typedef typename BaseType :: CoordinateType CoordinateType;

      //! type of key to identify a selected quadrature
      typedef typename QuadratureRuleRegistryType :: QuadratureKeyType  QuadratureKeyType;

    protected:
      //enum { highest_order_cube    = CubeQuadratureRule<ct,dim>::highest_order };
      //enum { highest_order_simplex = SimplexQuadratureRule<ct,dim>::highest_order };

      //(highest_order_cube < highest_order_simplex) ?
      //            highest_order_cube : highest_order_simplex };

    protected:
      const GeometryType elementGeometry_;
      int order_;

    public:
      /** \brief constructor filling the list of points and weights
       *
       *  \param[in]  geometry  geometry type for which a quadrature is desired
       *  \param[in]  key       key to identify quadrature (i.e. order of quadrature provided by the user)
       *  \param[in]  id        unique identifier (provided by QuadratureProvider)
       */
      FempyQuadratureRulesFactory( const GeometryType &geometry,
                                   const QuadratureKeyType& key,
                                   const size_t id )
      : BaseType( id ),
        elementGeometry_( geometry )
      {

        // if quadrature for this key was registered, then use it
        if( QuadratureRuleRegistryType::quadratureExists( key ) )
        {
          const DuneQuadratureRuleType &rule =
            QuadratureRuleRegistryType::quadratureRule( key );

          // obtain quadrature order
          order_ = key.first();

          typedef typename DuneQuadratureRuleType :: iterator IteratorType;
          const IteratorType endit = rule.end();
          for( IteratorType it = rule.begin(); it != endit; ++it )
            addQuadraturePoint( (*it).position(), (*it).weight() );
        }
        else
        {
          // when no quadrature for a key is registered
          // use default dune-fem quadratures
          DefaultQuadrature quad( geometry, key.first(), id );

          // obtain quadrature order
          order_ = quad.order();

          const int nop = quad.nop();
          for( int i=0; i<nop; ++i )
          {
            addQuadraturePoint( quad.point( i ), quad.weight( i ) );
          }
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
        return QuadratureKeyType :: highest_order;
      }
    };


  } // namespace FemPy

} // namespace Dune

#endif // #ifndef DUNE_FEMPY_FEMPYQUADRATURES_HH
