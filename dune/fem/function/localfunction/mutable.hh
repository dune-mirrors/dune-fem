#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_MUTABLE_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_MUTABLE_HH

//-s system includes
#include <cassert>
#include <utility>

//- Dune includes
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{

  namespace Fem
  {

    template< class >
    struct DiscreteFunctionTraits;

    //**************************************************************************
    //
    //  --MutableLocalFunction
    //
    //**************************************************************************
    //! Manages the getting and deleting of local function pointers and
    //! acts like a local functions
    template < class DiscreteFunction >
    class MutableLocalFunction
    : public LocalFunction< typename DiscreteFunctionTraits< DiscreteFunction > :: DiscreteFunctionSpaceType ::  BasisFunctionSetType,
      typename DiscreteFunctionTraits< DiscreteFunction > :: LocalDofVectorType >
    {
      typedef MutableLocalFunction< DiscreteFunction > ThisType;
      typedef LocalFunction< typename DiscreteFunctionTraits< DiscreteFunction > :: DiscreteFunctionSpaceType ::  BasisFunctionSetType,
      typename DiscreteFunctionTraits< DiscreteFunction > :: LocalDofVectorType > BaseType;

    public:
      //! type of DiscreteFunction
      typedef DiscreteFunction DiscreteFunctionType;

      //! type of the entity, the local function lives on is given by the space
      typedef typename BaseType::EntityType EntityType;

      //! type of local Dof vector object
      typedef typename BaseType::LocalDofVectorType LocalDofVectorType;

      //! type of BasisFunctionSet
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;

      //! Constructor creating empty local function from given discrete function
      explicit MutableLocalFunction ( DiscreteFunctionType &discreteFunction )
      : BaseType( LocalDofVectorType( discreteFunction.localDofVectorAllocator() ) ),
        discreteFunction_( &discreteFunction )
      {}

      //! Constructor creating empty local function from given discrete function
      explicit MutableLocalFunction ( const DiscreteFunctionType &discreteFunction )
      : BaseType( LocalDofVectorType( discreteFunction.localDofVectorAllocator() ) ),
        discreteFunction_( &const_cast<DiscreteFunctionType &>( discreteFunction ) )
      {}

      //! Constructor creating local function from given discrete function and entity, not empty
      explicit MutableLocalFunction ( DiscreteFunctionType &discreteFunction, const EntityType &entity )
      : BaseType( discreteFunction.space().basisFunctionSet( entity ), LocalDofVectorType( discreteFunction.localDofVectorAllocator() ) ),
        discreteFunction_( &discreteFunction )
      {
        discreteFunction.getLocalDofReferences( entity, localDofVector() );
      }

      //! Constructor creating local function from given discrete function and entity, not empty
      explicit MutableLocalFunction ( const DiscreteFunctionType &dFunction, const EntityType &entity )
      : BaseType( dFunction.space().basisFunctionSet( entity ), LocalDofVectorType( dFunction.localDofVectorAllocator() ) ),
        discreteFunction_( &const_cast<DiscreteFunctionType &>( dFunction ) )
      {
        discreteFunction().getLocalDofReferences( entity, localDofVector() );
      }

      //! copy constructor
      MutableLocalFunction ( const ThisType &other )
      : BaseType( static_cast< const BaseType& > ( other ) ), discreteFunction_( other.discreteFunction_ )
      {}

      //! move constructor
      MutableLocalFunction ( ThisType &&other )
      : BaseType( static_cast< BaseType&& > ( other ) ), discreteFunction_( other.discreteFunction_ )
      {}

      ThisType& operator= ( const ThisType& ) = delete;
      ThisType& operator= ( ThisType&& ) = delete;

      using BaseType::localDofVector;

      void init ( const EntityType &entity )
      {
        BaseType::init( discreteFunction().space().basisFunctionSet( entity ) );
        discreteFunction().getLocalDofReferences( entity, localDofVector() );
      }

      void bind ( const EntityType &entity )
      {
        init(entity);
      }

      void unbind() { BaseType::unbind(); }

      const DiscreteFunctionType &discreteFunction () const
      {
        return *discreteFunction_;
      }
      DiscreteFunctionType &discreteFunction ()
      {
        return *discreteFunction_;
      }

    protected:
      DiscreteFunctionType *discreteFunction_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_MUTABLE_HH
