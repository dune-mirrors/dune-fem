#ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH
#define DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH

#include <utility>
#include <dune/common/dynvector.hh>
#include <dune/fem/function/localfunction/mutable.hh>
#include <dune/fem/function/localfunction/localfunction.hh>

namespace Dune
{

  namespace Fem
  {

    // external forward declerations
    // -----------------------------

    template<class>
    struct DiscreteFunctionTraits;

    template<class>
    class MutableLocalFunction;


    // internal forward declerations
    // -----------------------------

    template< class DiscreteFunction >
    class ConstLocalFunction;


    template < class BasisFunctionSet, class LocalDofVector >
    class BasicConstLocalFunction
    : public LocalFunction< BasisFunctionSet, LocalDofVector >
    {
      typedef BasicConstLocalFunction< BasisFunctionSet, LocalDofVector >  ThisType;
      typedef LocalFunction< BasisFunctionSet, LocalDofVector > BaseType;

    public:
      //! type of Dof
      typedef typename BaseType::DofType DofType;

      //! type of Entity
      typedef typename BaseType :: EntityType EntityType;

      //! type of BasisFunctionSet
      typedef typename BaseType :: BasisFunctionSetType BasisFunctionSetType;

      //! type of LocalDofVector
      typedef typename BaseType :: LocalDofVectorType LocalDofVectorType;

      //! default ctor
      BasicConstLocalFunction () {}

      explicit BasicConstLocalFunction ( const BasisFunctionSetType & basisFunctionSet ) : BaseType( basisFunctionSet ) {}
     
      explicit BasicConstLocalFunction ( const LocalDofVectorType &localDofVector ) : BaseType( localDofVector ) {}

      BasicConstLocalFunction ( const BasisFunctionSetType &basisFunctionSet, const LocalDofVectorType &localDofVector )
      : BaseType( basisFunctionSet, localDofVector )
      {}

      explicit BasicConstLocalFunction ( LocalDofVectorType &&localDofVector ) : BaseType( localDofVector ) {}

      BasicConstLocalFunction ( const BasisFunctionSetType &basisFunctionSet, LocalDofVectorType &&localDofVector )
      : BaseType( basisFunctionSet, localDofVector )
      {}

      BasicConstLocalFunction ( const BaseType &other ) : BaseType( other.basisFunctionSet(), other.localDofVector() ) {}

      BasicConstLocalFunction ( const ThisType & ) = default;
      BasicConstLocalFunction ( ThisType && ) = default;

      using BaseType::operator[];
      using BaseType::localDofVector;

   protected:
      DofType &operator[] ( int num ) 
      {
        return static_cast< BaseType &>( *this )[ num ];
      }

      using BaseType::clear;
      using BaseType::assign; 
      using BaseType::operator +=;
      using BaseType::operator -=;
      using BaseType::axpy;
    };

    /** \ingroup LocalFunction
        \class ConstLocalFunction
        \brief A constant local function carrying values for one entity
      
        A ConstLocalFunction is a LocalFunction which is basically doing the same as the 
        LocalFunction of a discrete function. The difference is that the local dofs 
        are not kept as references but are copied to a local storage. 
        Therefore, this is a const local function and any modification of dofs is not
        allowed. 
      
        \note Local DoF numbers correspond directly to array indices. Hence it
        may be more cache efficient to generate a ConstLocalFunction when only a 
        const access to the local function is needed. 
      
        \param DiscreteFunction type of the discrete function, the
                                local function shall belong to
     */
    template< class DiscreteFunction >
    class ConstLocalFunction
    : public BasicConstLocalFunction<
      typename DiscreteFunctionTraits<
      typename remove_const< DiscreteFunction > :: type >::DiscreteFunctionSpaceType::BasisFunctionSetType,
      Dune::DynamicVector< typename DiscreteFunctionTraits< typename remove_const< DiscreteFunction > :: type >::DofType,
        typename DiscreteFunctionTraits< typename remove_const< DiscreteFunction > :: type >::LocalDofVectorAllocatorType
      :: template rebind< typename DiscreteFunctionTraits< typename remove_const< DiscreteFunction >::type > ::DofType > ::other > >
    {
      typedef ConstLocalFunction< DiscreteFunction > ThisType;
      typedef BasicConstLocalFunction< typename DiscreteFunctionTraits< typename remove_const< DiscreteFunction > :: type >::DiscreteFunctionSpaceType::BasisFunctionSetType,
              Dune::DynamicVector< typename DiscreteFunctionTraits< typename remove_const< DiscreteFunction > :: type >::DofType,
              typename DiscreteFunctionTraits< typename remove_const< DiscreteFunction > :: type > :: LocalDofVectorAllocatorType
              :: template rebind< typename DiscreteFunctionTraits< typename remove_const< DiscreteFunction >::type >::DofType >::other  > >
          BaseType;

    public:
      typedef typename remove_const< DiscreteFunction >::type DiscreteFunctionType;
      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

      typedef typename BaseType::DofType DofType;
      typedef typename BaseType::EntityType EntityType;
      typedef typename BaseType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename BaseType::LocalDofVectorType LocalDofVectorType;

      /** \brief constructor creating a local function without binding it to an 
                 entity
        
          Creates the local function without initializing the fields depending on
          the current entity.
        
          \note Before using the local function it must be initilized by
          \code
          localFunction.init( entity );
          \endcode
        
          \param[in] df discrete function the local function shall belong to
       */
      explicit ConstLocalFunction ( const DiscreteFunctionType &df )
      : BaseType( LocalDofVectorType( df.localDofVectorAllocator() ) ),
        discreteFunction_( &df )
      {
      }

      //! cast a MutableLocalFunction into this one !!! expensive !!!
      ConstLocalFunction ( const typename DiscreteFunctionType::LocalFunctionType &localFunction )
      : BaseType( localFunction.basisFunctionSet(), LocalDofVectorType( localFunction.size(), localFunction.discreteFunction().localDofVectorAllocator() ) ),
        discreteFunction_( &localFunction.discreteFunction() )
      {
        std::copy( localFunction.localDofVector().begin(), localFunction.localDofVector().end(), localDofVector().begin() );
      }
      
      /** \brief constructor creating a local function and binding it to an
                 entity
        
          Creates the local function and initilizes the fields depending on the
          current entity. It is not necessary, though allowed, to call init
          before using the discrete function.
        
          \note The degrees of freedom are not initialized by this function.
          
          \param[in] df      discrete function the local function shall
                             belong to
          \param[in] entity  entity for initialize the local function to
       */
      ConstLocalFunction ( const DiscreteFunctionType &df, const EntityType &entity )
      : BaseType( df.space().basisFunctionSet( entity ), LocalDofVectorType( df.localDofVectorAllocator() )  ),
        discreteFunction_( &df )
      {
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }

      ConstLocalFunction ( const ThisType &other ) = default;
      ConstLocalFunction ( ThisType &&other ) = default;

      using BaseType::localDofVector;

      /** \copydoc Dune::Fem::LocalFunction :: init */
      void init ( const EntityType &entity )
      {
        BaseType::init( discreteFunction().space().basisFunctionSet( entity ) );
        discreteFunction().getLocalDofs( entity, localDofVector() );
      }

      const DiscreteFunctionType& discreteFunction() const { return *discreteFunction_; }

    protected:
      const DiscreteFunctionType* discreteFunction_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_LOCALFUNCTION_CONST_HH
