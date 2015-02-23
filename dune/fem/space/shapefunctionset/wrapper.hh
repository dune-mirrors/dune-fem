#ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_WRAPPER_HH
#define DUNE_FEM_SPACE_SHAPEFUNCTIONSET_WRAPPER_HH

// C++ includes
#include <cstddef>

/*
  @file
  @author Christoph Gersbacher
  @brief Identical wrapper for implementations of Dune::Fem::ShapeFunctionSet
*/


namespace Dune
{

  namespace Fem
  {

    // IdShapeFunctionSet
    // ------------------

    template< class ShapeFunctionSet >
    class IdShapeFunctionSet
    {
    public:
      typedef typename ShapeFunctionSet::FunctionSpaceType FunctionSpaceType;

      typedef typename ShapeFunctionSet::DomainType DomainType;
      typedef typename ShapeFunctionSet::RangeType RangeType;
      typedef typename ShapeFunctionSet::JacobianRangeType JacobianRangeType;
      typedef typename ShapeFunctionSet::HessianRangeType HessianRangeType;

      IdShapeFunctionSet ( const ShapeFunctionSet &shapeFunctionSet = ShapeFunctionSet() )
      : shapeFunctionSet_( shapeFunctionSet )
      {}

      int order () const
      {
        return shapeFunctionSet_.order();
      }

      std::size_t size () const
      {
        return shapeFunctionSet_.size();
      }

      template< class Point, class Functor >
      void evaluateEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.evaluateEach( x, functor );
      }

      template< class Point, class Functor >
      void jacobianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.jacobianEach( x, functor );
      }

      template< class Point, class Functor >
      void hessianEach ( const Point &x, Functor functor ) const
      {
        return shapeFunctionSet_.hessianEach( x, functor );
      }

    private:
      ShapeFunctionSet shapeFunctionSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_SHAPEFUNCTIONSET_WRAPPER_HH
