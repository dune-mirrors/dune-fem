#ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDRESHAPEFUNCTIONSETS_HH
#define DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDRESHAPEFUNCTIONSETS_HH

// C++ includes
#include <cstddef>

// dune-common includes
#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/common/nullptr.hh>

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

  // LegendreShapeFunctionSet
  // ------------------------

  template< class FunctionSpace, int polOrder >
  class LegendreShapeFunctionSet
  : public Dune::Fem::ShapeFunctionSet< FunctionSpace, LegendreShapeFunctionSet< FunctionSpace, polOrder > >
  {
    dune_static_assert( (FunctionSpace::dimRange == 1),
                        "FunctionSpace must be scalar (i.e., dimRange = 1)." );

    typedef LegendreShapeFunctionSet< FunctionSpace, polOrder > ThisType;
    typedef Dune::Fem::ShapeFunctionSet< FunctionSpace, ThisType > BaseType;

    static const int dimension = FunctionSpace::dimDomain;
    typedef Dune::Fem::LegendreDGBaseFunction< FunctionSpace, polOrder > ImplementationType;

  public:
    typedef typename BaseType::FunctionSpaceType FunctionSpaceType;
    typedef typename BaseType::DomainType DomainType;
    typedef typename BaseType::RangeType RangeType;
    typedef typename BaseType::JacobianRangeType JacobianRangeType;
    typedef typename BaseType::HessianRangeType HessianRangeType;

    static const int numShapeFunctions = Dune::Fem::NumLegendreBaseFunctions< polOrder, dimension >::numBaseFct;

    LegendreShapeFunctionSet ( const Dune::GeometryType &type )
    : type_( type ),
      shapeFunctions_( numShapeFunctions, nullptr )
    {
      for( int i = 0; i < numShapeFunctions; ++i )
        shapeFunctions_[ i ] = new ImplementationType( i );
    }
    
    ~LegendreShapeFunctionSet ()
    {
      for( int i = 0; i < numShapeFunctions; ++i )
      {
        if( shapeFunctions_[ i ] )
          delete shapeFunctions_[ i ];
      }
    }

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
      for( std::size_t i = 0; i < size(); ++i )
      {
        RangeType val;
        evaluate( i, coordinate( x ), val );
        functor( i, val );
      }
    }

    template< class Point, class Functor >
    void jacobianEach ( const Point &x, Functor functor ) const
    {
      for( std::size_t i = 0; i < size(); ++i )
      {
        JacobianRangeType jac;
        jacobian( i, coordinate( x ), jac );
        functor( i, jac );
      }
    }

    template< class Point, class Functor >
    void hessianEach ( const Point &x, Functor functor ) const
    {
      DUNE_THROW( Dune::NotImplemented, "Method hessianEach() not implemented yet" );
    }

  private:
    void evaluate( const int baseFunction, const DomainType &x, RangeType &value ) const
    {
      assert( baseFunction >= 0 && baseFunction < numShapeFunctions );
      Dune::FieldVector< int, 0 > diffVariable;
      shapeFunctions_[ baseFunction ].evaluate( diffVariable, x, value );
    }

    void jacobian ( const int baseFunction, const DomainType &x, JacobianRangeType &jacobian ) const
    {
      assert( baseFunction >= 0 && baseFunction < numShapeFunctions );
      Dune::FieldVector< int, 1 > diffVariable;
      int &i = diffVariable[ 0 ];
      for( i = 0; i < dimension; ++i )
      {
        RangeType tmp;
        shapeFunctions_[ baseFunction ].evaluate( diffVariable, x, tmp );
        jacobian[ 0 ][ i ] = tmp[ 0 ];
      }
    }

    Dune::GeometryType type_; 
    Dune::array< ImplementationType *, numShapeFunctions > shapeFunctions_;
  };

} // namespace AnisotropicDG 

#endif // #ifndef DUNE_FEM_SPACE_ANISOTROPICDGSPACE_LEGENDRESHAPEFUNCTIONSETS_HH
