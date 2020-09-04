#ifndef DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH
#define DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH

// C++ includes
#include <cstddef>
#include <vector>

#include <dune/fem/common/coordinate.hh>

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

      const ThisType *clone () const = 0;
    };



    // SimpleShapeFunctionSet
    // ----------------------

    template< class ShapeFunction >
    class SimpleShapeFunctionSet
    {
      typedef SimpleShapeFunctionSet< ShapeFunction > ThisType;

    public:
      typedef ShapeFunction ShapeFunctionType;

      typedef ThisType ScalarFunctionSpaceType;

      typedef typename ShapeFunction::FunctionSpaceType FunctionSpaceType;
      typedef typename FunctionSpaceType::DomainType DomainType;
      typedef typename FunctionSpaceType::RangeType RangeType;
      typedef typename FunctionSpaceType::JacobianRangeType JacobianRangeType;
      typedef typename FunctionSpaceType::HessianRangeType HessianRangeType;

      // this number is positive if the shape function set
      // is using Lagrange polynomials
      static const int lagrangePointId = -1;

      template< class Factory >
      explicit SimpleShapeFunctionSet ( const Factory &factory );

      SimpleShapeFunctionSet ( const ThisType &other );

      const ThisType &operator= ( const ThisType &other );

      ~SimpleShapeFunctionSet ();

      int order () const { return order_; }

      // Shape Function Set Interface Methods
      std::size_t size () const { return shapeFunctions_.size(); }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const;

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const;

    protected:
      std::vector< const ShapeFunctionType * > shapeFunctions_;
      int order_;
    };



    // Implementation of SimpleShapeFunctionSet
    // ----------------------------------------

    template< class ShapeFunction >
    template< class Factory >
    inline SimpleShapeFunctionSet< ShapeFunction >
      ::SimpleShapeFunctionSet ( const Factory &factory )
    {
      const std::size_t numShapeFunctions = factory.numShapeFunctions();
      shapeFunctions_.resize( numShapeFunctions );
      for( std::size_t i = 0; i < numShapeFunctions; ++i )
        shapeFunctions_[ i ] = factory.createShapeFunction( i );
      order_ = factory.order();
    }

    template< class ShapeFunction >
    inline SimpleShapeFunctionSet< ShapeFunction >
      ::SimpleShapeFunctionSet ( const ThisType &other )
    {
      *this = other;
    }

    template< class ShapeFunction >
    inline const typename SimpleShapeFunctionSet< ShapeFunction >::ThisType &
    SimpleShapeFunctionSet< ShapeFunction >::operator= ( const ThisType &other )
    {
      if( this == &other )
        return *this;

      for( std::size_t i = 0; i < size(); ++i )
        delete shapeFunctions_[ i ];

      const std::size_t numShapeFunctions = other.size();
      shapeFunctions_.resize( numShapeFunctions );
      for( std::size_t i = 0; i < numShapeFunctions; ++i )
        shapeFunctions_[ i ] = other.shapeFunctions_[ i ]->clone();

      order_ = other.order_;
      return *this;
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
        shapeFunctions_[ i ]->jacobian( coordinate( x ), jacobian );
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
        shapeFunctions_[ i ]->hessian( coordinate( x ), hessian );
        functor( i, hessian );
      }
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SHAPEFUNCTIONSET_SIMPLE_HH
