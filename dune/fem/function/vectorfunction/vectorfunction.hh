#ifndef DUNE_FEM_VECTORFUNCTION_HH
#define DUNE_FEM_VECTORFUNCTION_HH

#include <memory>
#include <string>

#include <dune/fem/common/referencevector.hh>
#include <dune/fem/common/stackallocator.hh>
#include <dune/fem/function/blockvectors/defaultblockvectors.hh>
#include <dune/fem/function/common/discretefunction.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/storage/envelope.hh>
#include <dune/fem/storage/vector.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template < class DiscreteFunctionSpace, class Vector >
    class VectorDiscreteFunction;

    // VectorDiscreteFunctionTraits
    // ----------------------------

    template< typename DiscreteFunctionSpace, typename Vector >
    struct DiscreteFunctionTraits< VectorDiscreteFunction< DiscreteFunctionSpace, Vector > >
      : public DefaultDiscreteFunctionTraits< DiscreteFunctionSpace,
                  SimpleBlockVector< Vector, DiscreteFunctionSpace::localBlockSize > >
    {
      typedef VectorDiscreteFunction< DiscreteFunctionSpace, Vector >  DiscreteFunctionType;
      typedef MutableLocalFunction< DiscreteFunctionType > LocalFunctionType;
    };


    // VectorDiscreteFunction
    // ----------------------

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

      VectorDiscreteFunction( const std::string &name,
                              const DiscreteFunctionSpaceType &space,
                              VectorType& vector )
        : BaseType( name, space ),
          vec_(),
          dofVector_( vector )
      {
      }

      VectorDiscreteFunction( const VectorDiscreteFunction& other )
        : BaseType( "copy of " + other.name(), other.space() ),
          vec_(),
          dofVector_( allocateDofVector( other.space() ) )
      {
        assign( other );
      }

      DofVectorType& dofVector() { return dofVector_; }
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
      // dof vector that stores referenc to vector
      DofVectorType dofVector_;
    };

    // Capabilibies
    // ------------

    namespace Capabilities
    {

      template< class DiscreteFunctionSpace, class DofVector >
      struct HasLeakPointer
        < Fem :: VectorDiscreteFunction< DiscreteFunctionSpace, DofVector > >
      : public HasLeakPointer< DofVector >
      {};

    }

  } // namespace Fem

} // namespace Dune

#include "managedvectorfunction.hh"

#endif // #ifndef DUNE_FEM_VECTORFUNCTION_HH
