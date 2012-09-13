#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDREDGBASISFUNCTIONS_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDREDGBASISFUNCTIONS_HH

// C++ includes
#include <cassert>

// dune-common includes
#include <dune/common/array.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/forloop.hh>

// dune-geometry includes
#include <dune/geometry/type.hh>

// dune-fem includes
#include <dune/fem/space/dgspace/legendredgbasefunctions.hh>
#include <dune/fem/space/shapefunctionset/shapefunctionset.hh>

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
*/


namespace AnisotropicDG
{

  // LegendreDGShapeFunctionSet
  // --------------------------

  template< class ScalarFunctionSpace, int polOrder >
  class LegendreDGShapeFunctionSet
  : public Dune::Fem::ShapeFunctionSet< ScalarFunctionSpace, LegendreDGShapeFunctionSet< ScalarFunctionSpace, polOrder > >
  {
    typedef LegendreDGShapeFunctionSet< ScalarFunctionSpace, polOrder > ThisType;
    typedef Dune::Fem::ShapeFunctionSet< ScalarFunctionSpace, ThisType > BaseType;

  public:
    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::HessianRangeType HessianRangeType;

    static const int dimension = FunctionSpaceType::dimDomain;
    static const int numShapeFunctions = Dune::Fem::NumLegendreBaseFunctions< polOrder, dimension >::numBaseFct;

    LegendreDGShapeFunctionSet ( const Dune::GeometryType &type )
    : type_( type )
    {}

    Dune::GeometryType type () const
    {
      return type_;
    }

    std::size_t size () const
    {
      return numShapeFunctions;
    }

    template< class Point, class Functor >
    void evaluateEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method evaluateEach() not implemented yet" );
    }

    template< class Point, class Functor >
    void jacobianEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method jacobianEach() not implemented yet" );
    }

    template< class Point, class Functor >
    void hessianEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method hessianEach() not implemented yet" );
    }

  private:
    Dune::GeometryType type_; 
  };

} // namespace AnisotropicDG 

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDREDGBASISFUNCTIONS_HH
