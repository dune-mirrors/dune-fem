#ifndef DUNE_FEM_VECTORFUNCTION_HH
#define DUNE_FEM_VECTORFUNCTION_HH

#include <memory>
#include <string>

#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/storage/envelope.hh>

namespace Dune
{
  namespace Fem
  {

    template < class DiscreteFunctionSpace, class Vector >
    class VectorDiscreteFunction;

#if HAVE_PETSC
    template <class DiscreteFunctionSpace>
    class PetscDiscreteFunction;
#endif


    template< typename DiscreteFunctionSpace, typename Vector >
    struct DiscreteFunctionTraits< VectorDiscreteFunction< DiscreteFunctionSpace, Vector > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace,
                  SimpleBlockVector< Vector, DiscreteFunctionSpace::localBlockSize > >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, Vector >  DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };



    template < class DiscreteFunctionSpace, class Vector >
    class VectorDiscreteFunction
    : public DiscreteFunctionDefault<
          VectorDiscreteFunction< DiscreteFunctionSpace, Vector > >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, Vector > ThisType;
      typedef DiscreteFunctionDefault< ThisType > BaseType;

    public:
      typedef Vector VectorType;
      typedef typename BaseType :: DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef typename BaseType :: DofVectorType              DofVectorType;
      typedef typename DofVectorType :: DofContainerType      DofContainerType;
      typedef typename BaseType :: DofType                    DofType;

      using BaseType::assign;

      /** \brief Constructor to use if the vector storing the dofs already exists
       *
       *  \param[in]  name         name of the discrete function
       *  \param[in]  space        space the discrete function lives in
       *  \param[in]  dofVector    reference to the dof vector
       */
      VectorDiscreteFunction( const std::string& name,
                              const DiscreteFunctionSpaceType& space,
                              VectorType& dofVector )
        : BaseType( name, space ),
          vec_(),
          dofVector_( dofVector )
      {}

      /** \brief Copy constructor */
      VectorDiscreteFunction( const ThisType& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          vec_(),
          dofVector_( allocateDofVector( other.space() ) )
      {
        assign( other );
      }

      /** \brief Move constructor */
      VectorDiscreteFunction( ThisType&& other )
        : BaseType( static_cast< BaseType && >( other ) ),
          vec_( std::move( other.vec_ ) ),
          dofVector_( other.dofVector_ )
      {}

      VectorDiscreteFunction () = delete;
      ThisType& operator= ( const ThisType& ) = delete;
      ThisType& operator= ( ThisType&& ) = delete;

#if HAVE_PETSC
      void assign( const PetscDiscreteFunction< DiscreteFunctionSpaceType >& g )
      {
        g.dofVector().copyTo( dofVector() );
      }
#endif

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      DofVectorType& dofVector() { return dofVector_; }

      /** \copydoc Dune::Fem::DiscreteFunctionInterface::dofVector() */
      const DofVectorType& dofVector() const { return dofVector_; }

    protected:
      // allocate managed dof storage
      VectorType& allocateDofVector ( const DiscreteFunctionSpaceType& space )
      {
        vec_.reset( new VectorType( space.size() ) );
        return *vec_;
      }

      // pointer to DofContainer
      std::unique_ptr< VectorType > vec_;
      // dof vector that stores reference to vector
      DofVectorType dofVector_;
    };

  } // namespace Fem
} // namespace Dune

#include "managedvectorfunction.hh"

#endif // #ifndef DUNE_FEM_VECTORFUNCTION_HH
