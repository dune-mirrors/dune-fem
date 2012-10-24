#ifndef DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTION_HH
#define DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTION_HH

// C++ includes
#include <cstdlib>

// dune-common includes
#include <dune/common/exceptions.hh>
#include <dune/common/static_assert.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/common/functionspace.hh>
#include <dune/fem/space/shapefunctionset/vectorial.hh>

// local includes
#include <dune/fem/space/lagrangespace/genericlagrangepoints.hh>
#include <dune/fem/space/lagrangespace/genericbasefunctions.hh>


namespace Dune
{

  namespace Fem
  {

    // LagrangeShapeFunctionInterface
    // ------------------------------

    /* \brief abstract base class for Lagrange shape functions
     *
     * \tparam  FunctionSpace  Scalar function space
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

    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    class LagrangeShapeFunction
    : public LagrangeShapeFunctionInterface< FunctionSpace >
    {
      typedef LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder > ThisType;
      typedef LagrangeShapeFunctionInterface< FunctionSpace > BaseType;

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

 

    // LagrangeScalarShapeFunctionSet
    // ------------------------------

    template< class FunctionSpace, int polOrder >
    class LagrangeScalarShapeFunctionSet
    {
      template< class Topology >
      struct Switch
      {
        static const unsigned int topologyId = Topology::id;
        
        typedef typename GeometryWrapper< topologyId, FunctionSpace::dimDomain >
          ::GenericGeometryType GenericGeometryType;

        typedef LagrangeShapeFunction< FunctionSpace, GenericGeometryType, polOrder > 
          ShapeFunctionType;
        typedef typename ShapeFunctionType::GenericBaseFunctionType GenericBaseFunctionType;
       
        static const int size = GenericLagrangePoint< GeometryType, polOrder >::numLagrangePoints;
      };

    public:
      explicit LagrangeScalarShapeFunctionSet ( const GeometryType &type );

      std::size_t size () const;

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const;
    };



    // LagrangeShapeFunctionSetTraits
    // ------------------------------

    template< class FunctionSpace, int order >
    struct LagrangeShapeFunctionSetTraits
    {
      typedef FunctionSpace FunctionSpaceType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename ToScalarFunctionSpace< FunctionSpaceType >::Type ScalarFunctionSpaceType;
      typedef LagrangeScalarShapeFunctionSet< ScalarFunctionSpaceType, order > ScalarShapeFunctionSetType;
    };



    // LagrangeShapeFunctionSet
    // ------------------------

    template< class FunctionSpace, int order >
    class LagrangeShapeFunctionSet
    : public VectorialShapeFunctionSet< typename LagrangeShapeFunctionSetTraits< FunctionSpace, order>::ScalarShapeFunctionSetType,
                                        typename LagrangeShapeFunctionSetTraits< FunctionSpace, order>::RangeType
                                      >
    {
      typedef LagrangeShapeFunctionSet< FunctionSpace, order > ThisType;
      typedef LagrangeShapeFunctionSetTraits< FunctionSpace, order > Traits;
      typedef typename Traits::ScalarShapeFunctionSetType ScalarShapeFunctionSetType;
      typedef VectorialShapeFunctionSet< ScalarShapeFunctionSetType, typename Traits::RangeType > BaseType;

    public:
      LagrangeShapeFunctionSet ( const GeometryType &type );
    };



    // Implementation of LagrangeShapeFunction
    // ---------------------------------------

    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::evaluate ( const DomainType &x, RangeType &value ) const
    {
      DUNE_THROW( NotImplemented, "Method evaluate() not implemented yet." );
    }

    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::jacobian ( const DomainType &x, JacobianRangeType &value ) const
    {
      DUNE_THROW( NotImplemented, "Method jacobian() not implemented yet." );
    }

    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline void LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::hessian ( const DomainType &x, HessianRangeType &value ) const
    {
      DUNE_THROW( NotImplemented, "Method hessian() not implemented yet." );
    }

    template< class FunctionSpace, class GeometryType, unsigned int polOrder >
    inline int LagrangeShapeFunction< FunctionSpace, GeometryType, polOrder >
      ::order () const
    {
      return polOrder;
    }



    // Implementation of LagrangeScalarShapeFunctionSet
    // ------------------------------------------------

#if 0
    template< class FunctionSpace, int polOrder >
    inline LagrangeScalarShapeFunctionSet< FunctionSpace, polOrder >
      ::LagrangeScalarShapeFunctionSet ( const GeometryType &type )
    {}
#endif

    template< class FunctionSpace, int polOrder >
    inline std::size_t LagrangeScalarShapeFunctionSet< FunctionSpace, polOrder >
      ::size () const
    {
      DUNE_THROW( NotImplemented, "Method size() not implemented yet." );
    }

    template< class FunctionSpace, int polOrder >
    template< class Point, class Functor >
    inline void LagrangeScalarShapeFunctionSet< FunctionSpace, polOrder >
      ::evaluateEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( NotImplemented, "Method evaluateEach() not implemented yet." );
    }

    template< class FunctionSpace, int polOrder >
    template< class Point, class Functor >
    inline void LagrangeScalarShapeFunctionSet< FunctionSpace, polOrder >
      ::jacobianEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( NotImplemented, "Method jacobianEach() not implemented yet." );
    }

    template< class FunctionSpace, int polOrder >
    template< class Point, class Functor >
    inline void LagrangeScalarShapeFunctionSet< FunctionSpace, polOrder >
      ::hessianEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( NotImplemented, "Method hessianEach() not implemented yet." );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_LAGRANGE_SHAPEFUNCTION_HH
