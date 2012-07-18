#ifndef DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH

#include <vector>

#include <dune/geometry/type.hh>

namespace Dune
{

  namespace Fem
  {

    // AbstractShapeFunction
    // ---------------------

    template< class FunctionSpace >
    class AbstractShapeFunction
    {
      typedef AbstractShapeFunction< FunctionSpace > ThisType;

    public:
      typedef FunctionSpace FunctionSpaceType;

      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      virtual ~AbstractShapeFunction () {}

      virtual void evaluate ( const DomainType &x, RangeType &value ) const = 0;

      virtual void jacobian ( const DomainType &x, JacobianRangeType &jacobian ) const = 0;

      virtual void hessian ( const DomainType &x, HessianRangeType &hessian ) const = 0;
    };



    // SimpleShapeFunctionSet
    // ----------------------

    template< class ShapeFunction >
    class SimpleShapeFunctionSet
    {
      typedef SimpleShapeFunctionSet< ShapeFunction > ThisType;

    public:
      typedef ShapeFunction ShapeFunctionType;
      
      typedef typename ShapeFunctionType::FunctionSpaceType FunctionSpaceType;

      template< class Factory >
      explicit SimpleShapeFunctionSet ( const GeometryType &geometryType, const Factory &factory );

      ~SimpleShapeFunctionSet ();


      // Shape Function Set Interface Methods

      GeometryType geometryType () const { return geometryType_; }
      
      std::size_t size () const { return shapeFunctions_.size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const;
     
    protected:
      GeometryType geometryType_;
      std::vector< const ShapeFunctionType * > shapeFunctions_;
    };



    // Implementation of SimpleShapeFunctionSet
    // ----------------------------------------

    template< class ShapeFunction >
    template< class Factory >
    inline SimpleShapeFunctionSet< ShapeFunction >
      ::SimpleShapeFunctionSet ( const GeometryType &geometryType, const Factory &factory )
    : geometryType_( geometryType )
    {
      const std::size_t numShapeFunctions = factory.numShapeFunctions();
      shapeFunctions_.resize( numShapeFunctions );
      for( std::size_t i = 0; i < numShapeFunctions; ++i )
        shapeFunctions_[ i ] = factory.createShapeFunction( i );
    }


    template< class ShapeFunction >
    inline SimpleShapeFunctionSet< ShapeFunction >::~SimpleShapeFunctionSet ()
    {
      for( std::size_t i = 0; i < size(); ++i )
        delete shapeFunctions_[ i ];
    }


    template< class ShapeFunction >
    template< class Point, class Functor >
    inline void SimpleShapeFunctionSet< ShapeFunction >
      ::evaluateEach ( const Point &x, Functor functor ) const
    {
      for( std::size_t i = 0; i < size(); ++i )
      {
        typename ShapeFunctionType::RangeType value;
        shapeFunctions_[ i ]->evaluate( coordinate( x ), value );
        functor( i, value );
      }
    }
    

    template< class ShapeFunction >
    template< class Point, class Functor >
    inline void SimpleShapeFunctionSet< ShapeFunction >
      ::jacobianEach ( const Point &x, Functor functor ) const
    {
      for( std::size_t i = 0; i < size(); ++i )
      {
        typename ShapeFunctionType::JacobianRangeType jacobian;
        shapeFunctions_[ i ]->evaluate( coordinate( x ), jacobian );
        functor( i, jacobian );
      }
    }


    template< class ShapeFunction >
    template< class Point, class Functor >
    inline void SimpleShapeFunctionSet< ShapeFunction >
      ::hessianEach ( const Point &x, Functor functor ) const
    {
      for( std::size_t i = 0; i < size(); ++i )
      {
        typename ShapeFunctionType::HessianRangeType hessian;
        shapeFunctions_[ i ]->evaluate( coordinate( x ), hessian );
        functor( i, hessian );
      }
    }

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH
