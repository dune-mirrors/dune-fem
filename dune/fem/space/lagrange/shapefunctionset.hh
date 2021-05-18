#ifndef DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTIONSET_HH

#include <cassert>
#include <cstdlib>

#include <dune/common/fvector.hh>
#include <dune/fem/common/forloop.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>

#include "genericbasefunctions.hh"
#include "genericlagrangepoints.hh"

/*
  @file
  @brief Shape function set for Lagrange space
  @author Christoph Gersbacher
*/

namespace Dune
{

  namespace Fem
  {

    // LagrangeShapeFunctionInterface
    // ------------------------------

    /**
     * \brief abstract base class for Lagrange shape functions
     *
     * \tparam  FunctionSpace  scalar function space
     */
    template< class FunctionSpace >
    class LagrangeShapeFunctionInterface
    {
      typedef LagrangeShapeFunctionInterface< FunctionSpace > ThisType;

      static_assert( (FunctionSpace::dimRange == 1), "FunctionSpace must be scalar." );

    public:
      typedef FunctionSpace FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      virtual ~LagrangeShapeFunctionInterface () {}

      virtual void evaluate ( const DomainType &x, RangeType &value ) const = 0;
      virtual void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const = 0;
      virtual void hessian ( const DomainType &x, HessianRangeType &hessian ) const = 0;

      virtual int order () const = 0;

      virtual const ThisType *clone () const = 0;
    };



    // LagrangeShapeFunction
    // ---------------------

    /**
     * \brief implementation of Lagrange shape function
     *        using generic Lagrange shape functions
     *
     * \tparam  FunctionSpace  scalar function space
     * \tparam  GeometryType   generic geometry type wrapper
     * \tparam  polOrder       polynomial order
     */
    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    class LagrangeShapeFunction
    : public LagrangeShapeFunctionInterface< FunctionSpace >
    {
      typedef LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder > ThisType;
      typedef LagrangeShapeFunctionInterface< FunctionSpace > BaseType;

      static const int dimension = FunctionSpace::dimDomain;

    public:
      typedef GenericLagrangeBaseFunction< FunctionSpace, GeometryType, polOrder >
        GenericBaseFunctionType;

      typedef typename BaseType::DomainType DomainType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;
      typedef typename BaseType::HessianRangeType HessianRangeType;

      explicit LagrangeShapeFunction ( const GenericBaseFunctionType &genericShapeFunction )
      : genericShapeFunction_( genericShapeFunction )
      {}

      virtual void evaluate ( const DomainType &x, RangeType &value ) const;
      virtual void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const;
      virtual void hessian ( const DomainType &x, HessianRangeType &hessian ) const;

      virtual int order () const { return polOrder; }

      virtual const BaseType *clone () const { return new ThisType( *this ); }

    protected:
      GenericBaseFunctionType genericShapeFunction_;
    };



    // LagrangeShapeFunctionFactory
    // ----------------------------

    /**
     * \brief factory class
     *
     * \tparam  FunctionSpace  scalar function space
     * \tparam  polOrder       polynomial order
     */

    template< class FunctionSpace, int maxPolOrder >
    class LagrangeShapeFunctionFactory
    {
      static const int dimension = FunctionSpace::dimDomain;

    public:
      typedef LagrangeShapeFunctionInterface< FunctionSpace > ShapeFunctionType;

    private:
      template< GeometryType::Id geometryId>
      struct Switch;

    public:
      explicit LagrangeShapeFunctionFactory ( const Dune::GeometryType &type, const int order = maxPolOrder )
        : gt_(type),
          order_( order )
      {}

      int order () const;

      std::size_t numShapeFunctions () const;

      ShapeFunctionType *createShapeFunction( std::size_t i ) const;

    private:
      const Dune::GeometryType gt_;
      const int order_;
    };



    // LagrangeShapeFunctionSet
    // ------------------------

    /**
     * \brief Lagrange shape function set
     *
     * \tparam  FunctionSpace  function space
     * \tparam  polOrder       polynomial order
     */
    template< class FunctionSpace, int maxPolOrder >
    class LagrangeShapeFunctionSet
    : public SimpleShapeFunctionSet<
        typename LagrangeShapeFunctionFactory< FunctionSpace, maxPolOrder >::ShapeFunctionType
      >
    {
      static_assert( (FunctionSpace::dimRange == 1), "FunctionSpace must be scalar." );

      typedef LagrangeShapeFunctionFactory< FunctionSpace, maxPolOrder > ShapeFunctionFactoryType;
      typedef SimpleShapeFunctionSet< typename ShapeFunctionFactoryType::ShapeFunctionType > BaseType;

    public:
      LagrangeShapeFunctionSet ( const Dune::GeometryType &type, const int order = maxPolOrder )
      : BaseType( ShapeFunctionFactoryType( type, order ) )
      {
      }
    };



    // Implementation of LagrangeShapeFunction
    // ---------------------------------------

    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::evaluate ( const DomainType &x, RangeType &value ) const
    {
      FieldVector< int, 0 > diffVariable;
      genericShapeFunction_.evaluate( diffVariable, x, value );
    }


    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const
    {
      FieldVector< int, 1 > diffVariable;
      RangeType tmp;

      int &i = diffVariable[ 0 ];
      for( i = 0; i < dimension; ++i )
      {
        genericShapeFunction_.evaluate( diffVariable, x, tmp );
        jacobian[ 0 ][ i ] = tmp[ 0 ];
      }
    }


    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::hessian ( const DomainType &x, HessianRangeType &hessian ) const
    {
      FieldVector< int, 2 > diffVariable;
      RangeType tmp;

      int &i = diffVariable[ 0 ];
      for( i = 0; i < dimension; ++i )
      {
        // we use symmetrized evaluation of the hessian, since calling
        // evaluate is in general quite expensive
        int &j = diffVariable[ 1 ];
        for( j = 0; j < i; ++j )
        {
          genericShapeFunction_.evaluate( diffVariable, x, tmp );
          hessian[ 0 ][ i ][ j ] = hessian[ 0 ][ j ][ i ] = tmp[ 0 ];
        }

        assert( j == i );
        genericShapeFunction_.evaluate( diffVariable, x, tmp );
        hessian[ 0 ][ i ][ i ] = tmp[ 0 ];
      }
    }



    // LagrangeShapeFunctionFactory::Switch
    // ------------------------------------

    template< class FunctionSpace, int maxPolOrder >
    template< GeometryType::Id geometryId>
    struct LagrangeShapeFunctionFactory< FunctionSpace, maxPolOrder >::Switch
    {
      // get generic geometry type
      static constexpr GeometryType gt = geometryId;
      static const unsigned int topologyId = gt.id();
      typedef typename GeometryWrapper< topologyId, dimension >
        ::ImplType ImplType;

      template <int polOrd>
      struct CheckOrder
      {
        // type of scalar shape function for current polynomial order
        typedef LagrangeShapeFunction< FunctionSpace, ImplType, polOrd >   ShapeFunctionImpl;
        typedef typename ShapeFunctionImpl::GenericBaseFunctionType GenericBaseFunctionType;

        static void apply ( const int order, std::size_t &size )
        {
          if( order == polOrd )
          {
            size = GenericLagrangePoint< ImplType, polOrd >::numLagrangePoints;
          }
        }

        static void apply ( const std::size_t &i, const int order, ShapeFunctionType *&shapeFunction )
        {
          if( order == polOrd )
          {
            shapeFunction = new ShapeFunctionImpl( GenericBaseFunctionType( i ) );
          }
        }
      };

      static void apply ( const int order, std::size_t &size )
      {
        // loop over all possible polynomial order to find correct size
        Dune::Fem::ForLoop< CheckOrder, 0, maxPolOrder >::apply( order, size );
      }

      static void apply ( const std::size_t &i, const int order, ShapeFunctionType *&shapeFunction )
      {
        // loop over all possible polynomial order to create correct shape function
        Dune::Fem::ForLoop< CheckOrder, 0, maxPolOrder >::apply( i, order, shapeFunction );
      }
    };



    // Implementation of LagrangeShapeFunctionFactory
    // ----------------------------------------------

    template< class FunctionSpace, int maxPolOrder >
    const int LagrangeShapeFunctionFactory< FunctionSpace, maxPolOrder >::dimension;


    template< class FunctionSpace, int maxPolOrder >
    inline int LagrangeShapeFunctionFactory< FunctionSpace, maxPolOrder >
      ::order () const
    {
      return order_;
    }


    template< class FunctionSpace, int polOrder >
    inline std::size_t LagrangeShapeFunctionFactory< FunctionSpace, polOrder >
      ::numShapeFunctions () const
    {
      std::size_t numShapeFunctions( 0 );
      Dune::Impl::toGeometryTypeIdConstant<dimension>(gt_, [&](auto geometryTypeId) {
        Switch<decltype(geometryTypeId)::value>::apply(order_, numShapeFunctions);
      });
      return numShapeFunctions;
    }


    template< class FunctionSpace, int polOrder >
    inline typename LagrangeShapeFunctionFactory< FunctionSpace, polOrder >::ShapeFunctionType *
    LagrangeShapeFunctionFactory< FunctionSpace, polOrder >
      ::createShapeFunction( const std::size_t i ) const
    {
      ShapeFunctionType *shapeFunction( nullptr );
      Dune::Impl::toGeometryTypeIdConstant<dimension>(gt_, [&](auto geometryTypeId) {
        Switch<decltype(geometryTypeId)::value>::apply(i, order_,shapeFunction);
      });
      assert( shapeFunction );
      return shapeFunction;
    }



    // Extern Template Instantiations
    // ------------------------------

    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 1, 1 >, 1 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 2, 1 >, 1 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 3, 1 >, 1 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 1, 1 >, 1 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 2, 1 >, 1 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 3, 1 >, 1 >;

    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 1, 1 >, 2 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 2, 1 >, 2 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 3, 1 >, 2 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 1, 1 >, 2 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 2, 1 >, 2 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 3, 1 >, 2 >;

    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 1, 1 >, 3 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 2, 1 >, 3 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< double, double, 3, 1 >, 3 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 1, 1 >, 3 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 2, 1 >, 3 >;
    extern template class LagrangeShapeFunctionFactory< FunctionSpace< float, float, 3, 1 >, 3 >;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTIONSET_HH
