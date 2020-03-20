#ifndef DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_SCALARPRODUCTS_HH
#define DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_SCALARPRODUCTS_HH

#include <utility>

#include <dune/fem/common/utility.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/tuplediscretefunction/discretefunction.hh>


namespace Dune
{

  namespace Fem
  {

    // ParallelScalarProduct
    // ---------------------

    template< class ... DiscreteFunctions >
    class ParallelScalarProduct< TupleDiscreteFunction< DiscreteFunctions ... > >
    {
      typedef ParallelScalarProduct< TupleDiscreteFunction< DiscreteFunctions ... > > ThisType;

      typedef std::tuple< ParallelScalarProduct< DiscreteFunctions > ... > ParallelScalarProductTuple;

      typedef decltype ( std::index_sequence_for< DiscreteFunctions ... >() ) Sequence;

    public:
      typedef TupleDiscreteFunction< DiscreteFunctions ... > DiscreteFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;
      typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;

      ParallelScalarProduct ( const DiscreteFunctionSpaceType &space )
        : tuple_( createTuple( space, Sequence() ) )
      {}

      const DiscreteFunctionSpaceType &space () const { return space_; }

      RangeFieldType scalarProductDofs ( const DiscreteFunctionType &x, const DiscreteFunctionType &y ) const
      {
        return scalarProductDofs( x, y, Sequence() );
      }

      template< class OtherDiscreteFunction >
      RangeFieldType scalarProductDofs ( const DiscreteFunctionType &x, const OtherDiscreteFunction &y ) const
      {
        DUNE_THROW( NotImplemented, "Method scalarProductDofs ( DofVectorType, OtherDofVector ) not implemented" );
        return RangeFieldType( 0 );
      }

    protected:
      template< std::size_t ... I >
      RangeFieldType scalarProductDofs ( const DiscreteFunctionType &x, const DiscreteFunctionType &y, std::index_sequence< I ... > ) const
      {
        return Std::sum( std::get< I >( tuple_ ).scalarProductDofs( x.template subDiscreteFunction< I >(), y.template subDiscreteFunction< I >() ) ... );
      }

      template< std::size_t ... I >
      static ParallelScalarProductTuple createTuple ( const DiscreteFunctionSpaceType &space, std::index_sequence< I ... > )
      {
        return std::make_tuple( std::tuple_element< I, ParallelScalarProductTuple >::type( space.template subDiscreteFunctionSpace< I >() ) ... );
      }

      const DiscreteFunctionSpaceType &space_;
      ParallelScalarProductTuple tuple_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_SCALARPRODUCTS_HH
