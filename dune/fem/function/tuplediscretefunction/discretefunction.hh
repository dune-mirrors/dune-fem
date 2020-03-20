#ifndef DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_DISCRETEFUNCTION_HH
#define DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_DISCRETEFUNCTION_HH

#include <string>
#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/common/scalarproducts.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/function/tuplediscretefunction/functor.hh>

#include <dune/fem/function/tuplediscretefunction/dofvector.hh>
#include <dune/fem/space/combinedspace.hh>

namespace Dune
{

  namespace Fem
  {

    //! forward declaration
    template< class ... DiscreteFunctions >
    class TupleDiscreteFunction;

    // DiscreteFunctionTraits
    // ----------------------

    template< class ... DiscreteFunctions >
    struct DiscreteFunctionTraits< TupleDiscreteFunction< DiscreteFunctions ... > >
      : public DefaultDiscreteFunctionTraits<
          TupleDiscreteFunctionSpace< typename DiscreteFunctions::DiscreteFunctionSpaceType ... >,
          TupleDofVector< typename DiscreteFunctions::DofVectorType ... >
          >
    {
      typedef TupleDiscreteFunction< DiscreteFunctions ... > DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    // TupleDiscreteFunction
    // ---------------------

    template< class ... DiscreteFunctions >
    class TupleDiscreteFunction
      : public DiscreteFunctionDefault< TupleDiscreteFunction< DiscreteFunctions ... > >,
        public std::tuple< DiscreteFunctions ... >
    {
      typedef TupleDiscreteFunction< DiscreteFunctions ... > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

      typedef ParallelScalarProduct< ThisType > ScalarProductType;

      typedef std::tuple< DiscreteFunctions ... > DiscreteFunctionTuple;

      static_assert( sizeof ... ( DiscreteFunctions ) > 0, "TupleDiscreteFunction needs at least one DiscreteFunction." );

    public:
      typedef decltype ( std::index_sequence_for< DiscreteFunctions ... >() ) Sequence;

      typedef TupleDofVector< typename DiscreteFunctions::DofVectorType ... > DofVectorType;

      using BaseType::space;

      //! type for the discrete function space this function lives in
      typedef TupleDiscreteFunctionSpace< typename DiscreteFunctions::DiscreteFunctionSpaceType ... > DiscreteFunctionSpaceType;

      //! helper struct to get the type of the i-th sub function
      template< int i >
      struct SubDiscreteFunction
      {
        typedef typename std::tuple_element< i, DiscreteFunctionTuple >::type Type;
      };

      /** \brief Constructor to use if the vector storing the dofs (which is a block vector) already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  dfSpace      space the discrete function lives in
       *  \param[in]  blockVector  reference to the blockVector
       */
      TupleDiscreteFunction ( const std::string &name,
                              const DiscreteFunctionSpaceType &dfSpace,
                              DofVectorType &dofVector )
        : TupleDiscreteFunction ( name, dfSpace, dofVector, Sequence() ) {}

      /** \brief Constructor to use if the vector storing the dofs does not exist yet
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  dfSpace      space the discrete function lives in
       */
      TupleDiscreteFunction ( const std::string &name,
                              const DiscreteFunctionSpaceType &dfSpace )
        : TupleDiscreteFunction ( name, dfSpace, Sequence() ) {}

      // copy constructor
      TupleDiscreteFunction ( const ThisType &other )
        : TupleDiscreteFunction ( "copy of "+other.name(), other.space(), Sequence() )
      {
        dofVector_ = other.dofVector();
      }

      // move constructor
      TupleDiscreteFunction ( ThisType&& other )
        : BaseType ( static_cast< BaseType && >( other ) ),
          DiscreteFunctionTuple( std::move( other ) ),
          dofVector_( std::move( other.dofVector_ ) )
      {}

      TupleDiscreteFunction () = delete;
      ThisType &operator= ( const ThisType & ) = delete;

      DofVectorType &dofVector () { return dofVector_; }
      const DofVectorType &dofVector () const { return dofVector_; }

      template< int i >
      typename SubDiscreteFunction< i >::Type& subDiscreteFunction()
      {
        return std::get< i >( *this );
      }

      template< int i >
      const typename SubDiscreteFunction< i >::Type& subDiscreteFunction() const
      {
        return std::get< i >( *this );
      }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::enableDofCompression() */
      void enableDofCompression ()
      {
        Hybrid::forEach( Sequence{}, [ & ]( auto i ){  std::get< i >( *this ).enableDofCompression(); } );
      }

    protected:

      template< std::size_t ... I >
      TupleDiscreteFunction ( const std::string &name, const DiscreteFunctionSpaceType &space, std::index_sequence< I ... > )
        : BaseType( name, space ),
          DiscreteFunctionTuple(
            typename SubDiscreteFunction< I >::Type(
              name + "_comp_" + std::to_string( I ), space.template subDiscreteFunctionSpace< I >()
              ) ... ),
          dofVector_( std::get< I >( *this ).dofVector() ... )
      {}

      template< std::size_t ... I >
      TupleDiscreteFunction ( const std::string &name, const DiscreteFunctionSpaceType &space,
                              DofVectorType &dofVector, std::index_sequence< I ... > )
        : BaseType( name, space ),
          DiscreteFunctionTuple(
            typename SubDiscreteFunction< I >::Type(
              name + "_comp_" + std::to_string( I ), space.template subDiscreteFunctionSpace< I >(),
              std::get< I >( dofVector ) ) ... ),
          dofVector_( dofVector )
      {}

      DofVectorType dofVector_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_TUPLEDISCRETEFUNCTION_DISCRETEFUNCTION_HH
