#ifndef DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTIONSET_HH

// C++ includes
#include <cassert>
#include <cstdlib>

// dune-common includes
#include <dune/common/fvector.hh>
#include <dune/common/nullptr.hh>
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/simple.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include "genericlagrangepoints.hh"
#include "genericbasefunctions.hh"

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
      dune_static_assert( (FunctionSpace::dimRange == 1), "FunctionSpace must be scalar." );

    public:
      typedef FunctionSpace FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      virtual void evaluate ( const DomainType &x, RangeType &value ) const = 0;
      virtual void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const = 0;
      virtual void hessian ( const DomainType &x, HessianRangeType &hessian ) const = 0;

      virtual int order () const = 0;
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

      explicit LagrangeShapeFunction ( const GenericBaseFunctionType &genericShapeFunction );

      virtual void evaluate ( const DomainType &x, RangeType &value ) const;
      virtual void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const;
      virtual void hessian ( const DomainType &x, HessianRangeType &hessian ) const;

      int order () const;
      
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

    template< class FunctionSpace, int polOrder >
    class LagrangeShapeFunctionFactory 
    {
      static const int dimension = FunctionSpace::dimDomain;

    public:
      typedef LagrangeShapeFunctionInterface< FunctionSpace > ShapeFunctionType;

    private:
      template< class Topology >
      struct Switch;

    public:
      explicit LagrangeShapeFunctionFactory ( const Dune::GeometryType &type );

      std::size_t numShapeFunctions () const;

      ShapeFunctionType *createShapeFunction( std::size_t i ) const;

    private:
      unsigned int topologyId_;
    };



    // LagrangeShapeFunctionSetTraits
    // ------------------------------

    template< class FunctionSpace, int order >
    struct LagrangeShapeFunctionSetTraits
    {
      typedef FunctionSpace FunctionSpaceType;

      typedef typename FunctionSpaceType::RangeType RangeType;

      typedef typename ToScalarFunctionSpace< FunctionSpaceType >::Type ScalarFunctionSpaceType;
      typedef LagrangeShapeFunctionFactory< ScalarFunctionSpaceType, order > ScalarShapeFunctionFactoryType;
      typedef typename ScalarShapeFunctionFactoryType::ShapeFunctionType ScalarShapeFunctionType;
      typedef SimpleShapeFunctionSet< ScalarShapeFunctionType > ScalarShapeFunctionSetType;
    };



    // LagrangeShapeFunctionSet
    // ------------------------

    /**
     * \brief Lagrange shape function set
     *
     * \tparam  FunctionSpace  function space
     * \tparam  polOrder       polynomial order
     */
    template< class FunctionSpace, int order >
    class LagrangeShapeFunctionSet
    : public VectorialShapeFunctionSet< typename LagrangeShapeFunctionSetTraits< FunctionSpace, order>::ScalarShapeFunctionSetType,
                                        typename LagrangeShapeFunctionSetTraits< FunctionSpace, order>::RangeType
                                      >
    {
      typedef LagrangeShapeFunctionSetTraits< FunctionSpace, order> Traits;

      typedef VectorialShapeFunctionSet< typename Traits::ScalarShapeFunctionSetType,
                                         typename Traits::RangeType
                                       > BaseType;

      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      typedef typename Traits::ScalarShapeFunctionFactoryType ScalarShapeFunctionFactoryType;

    public:
      typedef typename BaseType::ScalarFunctionSpaceType ScalarFunctionSpaceType;

      LagrangeShapeFunctionSet ( const Dune::GeometryType &type );
    };



    // Implementation of LagrangeShapeFunction
    // ---------------------------------------

    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::evaluate ( const DomainType &x, RangeType &value ) const
    {
      FieldVector< int, 0 > diffVariable;
      return genericShapeFunction_.evaluate( diffVariable, x, value );
    }


    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::jacobian ( const DomainType &x, JacobianRangeType &value ) const
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
      ::hessian ( const DomainType &x, HessianRangeType &value ) const
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


    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline int LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::order () const
    {
      return polOrder;
    }



    // LagrangeShapeFunctionFactory::Switch
    // ------------------------------------

    template< class FunctionSpace, int polOrder >
    template< class Topology >
    struct LagrangeShapeFunctionFactory< FunctionSpace, polOrder >::Switch
    {
      // get generic geometry type
      static const unsigned int topologyId = Topology::id;
      typedef typename GeometryWrapper< topologyId, dimension >
        ::GenericGeometryType GenericGeometryType;

      // type of scalar shape function
      typedef LagrangeShapeFunction< FunctionSpace, GenericGeometryType, polOrder > 
        ShapeFunctionImpl;
      typedef typename ShapeFunctionImpl::GenericBaseFunctionType GenericBaseFunctionType;

      static void apply ( std::size_t &size )
      {
        size = GenericLagrangePoint< GenericGeometryType, polOrder >::numLagrangePoints;
      }

      static void apply ( const std::size_t &i, ShapeFunctionType *&shapeFunction )
      {
        shapeFunction = new ShapeFunctionImpl( GenericBaseFunctionType( i ) );
      }
    };



    // Implementation of LagrangeShapeFunctionFactory
    // ----------------------------------------------

    template< class FunctionSpace, int polOrder >
    inline LagrangeShapeFunctionFactory< FunctionSpace, polOrder >
      ::LagrangeShapeFunctionFactory ( const Dune::GeometryType &type )
    : topologyId_( type.id() )
    {}


    template< class FunctionSpace, int polOrder >
    inline std::size_t LagrangeShapeFunctionFactory< FunctionSpace, polOrder >
      ::numShapeFunctions () const
    {
      std::size_t numShapeFunctions;
      GenericGeometry::IfTopology< Switch, dimension >::apply( topologyId_, numShapeFunctions );
      return numShapeFunctions;
    }


    template< class FunctionSpace, int polOrder >
    inline typename LagrangeShapeFunctionFactory< FunctionSpace, polOrder >::ShapeFunctionType *
    LagrangeShapeFunctionFactory< FunctionSpace, polOrder >
      ::createShapeFunction( const std::size_t i ) const
    {
      ShapeFunctionType *shapeFunction( nullptr );
      GenericGeometry::IfTopology< Switch, dimension >::apply( topologyId_, i, shapeFunction );
      return shapeFunction;
    }



    // Implementation of LagrangeShapeFunctionSet
    // ------------------------------------------

    template< class FunctionSpace, int order >
    inline LagrangeShapeFunctionSet< FunctionSpace, order >
      ::LagrangeShapeFunctionSet( const Dune::GeometryType &type )
    : BaseType( ScalarShapeFunctionSetType( ScalarShapeFunctionFactoryType( type ) ) )
    {}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTIONSET_HH
