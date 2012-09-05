#ifndef DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH

#include <vector>

#include <dune/geometry/type.hh>

#include <dune/fem/space/shapefunctionset/shapefunctionset.hh>

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
    : public ShapeFunctionSet< typename ShapeFunction::FunctionSpaceType, SimpleShapeFunctionSet< ShapeFunction > >
    {
      typedef SimpleShapeFunctionSet< ShapeFunction > ThisType;
      typedef ShapeFunctionSet< typename ShapeFunction::FunctionSpaceType, 
                                SimpleShapeFunctionSet< ShapeFunction > > BaseType;

    public:
      typedef ShapeFunction ShapeFunctionType;
      
      typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
      typedef typename BaseType::DomainType DomainType;
      typedef typename BaseType::RangeType RangeType;
      typedef typename BaseType::JacobianRangeType JacobianRangeType;
      typedef typename BaseType::HessianRangeType HessianRangeType;

      template< class Factory >
      explicit SimpleShapeFunctionSet ( const GeometryType &type, const Factory &factory );

      ~SimpleShapeFunctionSet ();


      // Shape Function Set Interface Methods

      GeometryType type () const { return type_; }
      
      std::size_t size () const { return shapeFunctions_.size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor > 
      void hessianEach ( const Point &x, Functor functor ) const;
     
    protected:
      GeometryType type_;
      std::vector< const ShapeFunctionType * > shapeFunctions_;
    };



    // Implementation of SimpleShapeFunctionSet
    // ----------------------------------------

    template< class ShapeFunction >
    template< class Factory >
    inline SimpleShapeFunctionSet< ShapeFunction >
      ::SimpleShapeFunctionSet ( const GeometryType &type, const Factory &factory )
    : type_( type )
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
        RangeType value;
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
        JacobianRangeType jacobian;
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
        HessianRangeType hessian;
        shapeFunctions_[ i ]->evaluate( coordinate( x ), hessian );
        functor( i, hessian );
      }
    }

  } // namespace Fem

} // namespace Dune 

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH
