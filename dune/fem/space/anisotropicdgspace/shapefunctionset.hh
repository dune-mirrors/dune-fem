#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH

// C++ includes
#include <cassert>

// dune-common includes
#include <dune/common/exceptions.hh>

// dune-geometry includes
#include <dune/geometry/referenceelements.hh>

// dune-fem includes
#include <dune/fem/space/shapefunctionset/shapefunctionset.hh>

// local includes
#include "multiindexset.hh"

/**
  @file
  @author Christoph Gersbacher
  @brief Please doc me
*/


namespace AnisotropicDG 
{

  // ShapeFunctionSet
  // ----------------

  template< class FunctionSpace, int maxOrder >
  class ShapeFunctionSet
  : public Dune::Fem::ShapeFunctionSet< FunctionSpace, ShapeFunctionSet< FunctionSpace, maxOrder > >
  {
    typedef ShapeFunctionSet< FunctionSpace, maxOrder > ThisType;
    typedef Dune::Fem::ShapeFunctionSet< FunctionSpace, ThisType > BaseType;

  public:
    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::HessianRangeType HessianRangeType;

  private:
    static const int dimension = FunctionSpaceType::dimDomain;

  public:
    typedef typename MultiIndexSet< dimension, maxOrder >::MultiIndexType MultiIndexType;

    ShapeFunctionSet ( const Dune::GeometryType &type, const MultiIndexType &multiIndex )
    : multiIndex_( multiIndex )
    {
      assert( type == ThisType::type() );
    }

    static Dune::GeometryType type ()
    {
      return Dune::ReferenceElements< typename FunctionSpaceType::DomainFieldType, dimension >::cube(); 
    }

    std::size_t size () const
    {
      DUNE_THROW( Dune::NotImplemented, "Method size() not implemented yet" );
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

  protected:
    const MultiIndexType &multiIndex () const
    {
      return multiIndex_;
    }

  private:
    MultiIndexType multiIndex_;
  };

} // namespace AnisotropicDG 

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_SHAPEFUNCTIONSET_HH
