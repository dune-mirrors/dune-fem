#ifndef DUNE_FEM_FUNCTION_COMMON_LOCALCONTRIBUTION_HH
#define DUNE_FEM_FUNCTION_COMMON_LOCALCONTRIBUTION_HH

#include <algorithm>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/densevector.hh>
#include <dune/common/ftraits.hh>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/common/localcontribution.hh>
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/function/localfunction/temporary.hh>

namespace Dune
{

  namespace Fem
  {

    // External Forward Declarations
    // -----------------------------

    template< class >
    struct DiscreteFunctionTraits;

    class IsDiscreteFunction;



    namespace Assembly
    {

      namespace Global
      {

        // AddBase
        // -------

        template< class DiscreteFunction >
        struct AddBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        {
          typedef typename DiscreteFunction::DofType DofType;

          static void begin ( DiscreteFunction &df )
          {
            typedef typename DiscreteFunction::DiscreteFunctionSpaceType::LocalBlockIndices LocalBlockIndices;

            // clear auxiliary DoFs
            auto &dofVector = df.dofVector();
            for( const auto &auxiliaryDof : df.space().auxiliaryDofs() )
              Hybrid::forEach( LocalBlockIndices(), [ &dofVector, &auxiliaryDof ] ( auto &&j ) { dofVector[ auxiliaryDof ][ j ] = DofType( 0 ); } );
          }

          static void end ( DiscreteFunction &df ) { df.space().communicate( df, DFCommunicationOperation::Add() ); }
        };



        // SetBase
        // -------

        template< class DiscreteFunction >
        struct SetBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        {
          static void begin ( DiscreteFunction &df ) {}
          static void end ( DiscreteFunction &df ) { df.space().communicate( df, DFCommunicationOperation::Copy() ); }
        };

      } // namespace Global



      // AddBase
      // -------

      template< class DiscreteFunction >
      struct AddBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
      {
        typedef typename DiscreteFunction::DofType DofType;

        typedef Global::Add< DiscreteFunction > GlobalOperationType;

        template< class Entity, class LocalDofVector >
        void begin ( const Entity &entity, const DiscreteFunction &df, LocalDofVector &localDofVector ) const
        {
          std::fill( localDofVector.begin(), localDofVector.end(), DofType( 0 ) );
        }

        template< class Entity, class LocalDofVector >
        void end ( const Entity &entity, LocalDofVector &localDofVector, DiscreteFunction &df ) const
        {
          df.addLocalDofs( entity, localDofVector );
        }
      };



      // AddScaledBase
      // -------------

      template< class DiscreteFunction >
      struct AddScaledBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        : public AddBase< DiscreteFunction >
      {
        AddScaledBase ( typename DiscreteFunction::DofType factor ) : factor_( std::move( factor ) ) {}

        template< class Entity, class LocalDofVector >
        void end ( const Entity &entity, LocalDofVector &localDofVector, DiscreteFunction &df ) const
        {
          df.addScaledLocalDofs( entity, factor_, localDofVector );
        }

      private:
        typename DiscreteFunction::DofType factor_;
      };


      namespace detail
      {

        template< class DiscreteFunction, const bool getAndSet >
        struct SetAndSelectDFImpl
        {
          typedef typename DiscreteFunction::DofType DofType;

          typedef Global::Set< DiscreteFunction > GlobalOperationType;

          template< class Entity, class LocalDofVector >
          void begin ( const Entity &entity, const DiscreteFunction &df, LocalDofVector &localDofVector ) const
          {
            if constexpr ( getAndSet )
            {
              // obtain local dofs
              df.getLocalDofs ( entity, localDofVector );
            }
            else
            {
              // reset all dofs
              std::fill( localDofVector.begin(), localDofVector.end(), DofType( 0 ) );
            }
          }

          template< class Entity, class LocalDofVector >
          void end ( const Entity &entity, LocalDofVector &localDofVector, DiscreteFunction &df ) const
          {
            df.setLocalDofs( entity, localDofVector );
          }
        };
      }


      // SetBase
      // -------

      template< class DiscreteFunction >
      struct SetBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        : public detail::SetAndSelectDFImpl< DiscreteFunction, false >
      {};

      // SetSelectedBase
      // ---------------

      template< class DiscreteFunction >
      struct SetSelectedBase< DiscreteFunction, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
        : public detail::SetAndSelectDFImpl< DiscreteFunction, true >
      {};

    } // namespace Assembly
  }

  // consistency with Dune::DenseVector and DenseMatrix
  template< class DiscreteFunction, template< class > class AssemblyOperation >
  struct FieldTraits< Fem::LocalContribution< DiscreteFunction,  AssemblyOperation, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > > >
    : public FieldTraits< typename DiscreteFunction::DofType >
  {
  };


  namespace Fem
  {

    // LocalContribution for Discrete Functions
    // ----------------------------------------

    template< class DiscreteFunction, template< class > class AssemblyOperation >
    class LocalContribution< DiscreteFunction, AssemblyOperation, std::enable_if_t< std::is_base_of< Fem::IsDiscreteFunction, DiscreteFunction >::value > >
      : public TemporaryLocalFunction< typename DiscreteFunction::DiscreteFunctionSpaceType >
    {
      typedef typename DiscreteFunction::DiscreteFunctionSpaceType  DiscreteFunctionSpaceType;
      typedef LocalContribution< DiscreteFunction, AssemblyOperation >  ThisType;
      typedef TemporaryLocalFunction< DiscreteFunctionSpaceType >       BaseType;

    public:
      typedef DiscreteFunction DiscreteFunctionType;
      typedef AssemblyOperation< typename DiscreteFunctionTraits< DiscreteFunctionType >::DiscreteFunctionType > AssemblyOperationType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType::BasisFunctionSetType BasisFunctionSetType;
      typedef typename DiscreteFunctionType::DofType DofType;

      typedef typename DiscreteFunctionType::RangeType RangeType;
      typedef typename RangeType::field_type RangeFieldType;
      typedef typename DiscreteFunctionType::JacobianRangeType JacobianRangeType;

      typedef typename BaseType :: LocalDofVectorType  LocalDofVectorType;
      typedef typename LocalDofVectorType::size_type SizeType;

      typedef typename BasisFunctionSetType::EntityType EntityType;

      using BaseType::entity;
      using BaseType::localDofVector;
      using BaseType::axpy;

      template< class... Args >
      explicit LocalContribution ( DiscreteFunctionType &discreteFunction, Args &&... args )
        : BaseType( discreteFunction.space() ),
          discreteFunction_( discreteFunction ),
          assemblyOperation_( std::forward< Args >( args )... ),
          bound_( false )
      {
        discreteFunction.template beginAssemble< typename AssemblyOperationType::GlobalOperationType >();
      }

      LocalContribution ( const ThisType & ) = delete;
      LocalContribution ( ThisType && ) = delete;

      ~LocalContribution () { discreteFunction().template endAssemble< typename AssemblyOperationType::GlobalOperationType >(); }

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      const DiscreteFunctionType& discreteFunction () const { return discreteFunction_; }
      DiscreteFunctionType& discreteFunction () { return discreteFunction_; }

      void bind ( const EntityType &entity )
      {
        BaseType::bind( entity );
        bound_ = true;
        assemblyOperation_.begin( entity, discreteFunction(), localDofVector() );
      }

      void unbind ()
      {
        if (bound_)
        {
          // write back dofs to discrete function
          assemblyOperation_.end( entity(), localDofVector(), discreteFunction() );
          // unbind local contribution
          BaseType::unbind();
        }
      }

    protected:
      // LocalContribution is not a LocalFunction,
      // thus disable evaluate,jacobian and hessian methods
      using BaseType::evaluate;
      using BaseType::evaluateQuadrature;
      using BaseType::jacobian;
      using BaseType::hessian;

    protected:
      DiscreteFunctionType& discreteFunction_;
      AssemblyOperationType assemblyOperation_;
      bool bound_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_FUNCTION_COMMON_LOCALCONTRIBUTION_HH
