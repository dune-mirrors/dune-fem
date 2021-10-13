#ifndef DUNE_FEM_OPERATOR_COMMON_LOCALCONTRIBUTION_HH
#define DUNE_FEM_OPERATOR_COMMON_LOCALCONTRIBUTION_HH

#include <algorithm>
#include <type_traits>
#include <vector>
#include <utility>

#include <dune/common/densematrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/fem/common/utility.hh>
#include <dune/fem/function/common/localcontribution.hh>
#include <dune/fem/operator/common/operator.hh>
#include <dune/fem/operator/common/temporarylocalmatrix.hh>
#include <dune/fem/storage/rowreferencevector.hh>
#include <dune/fem/storage/subvector.hh>

namespace Dune
{

  namespace Fem
  {

    namespace Assembly
    {

      namespace Global
      {

        // AddBase
        // -------

        template< class AssembledOperator >
        struct AddBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
        {
          template <class DF, class RF>
          static void begin ( Dune::Fem::AssembledOperator< DF, RF > &op ) {}
          template <class DF, class RF>
          static void end ( Dune::Fem::AssembledOperator< DF, RF > &op ) { op.flushAssembly(); }
        };



        // SetBase
        // -------

        template< class AssembledOperator >
        struct SetBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
        {
          template <class DF, class RF>
          static void begin ( Dune::Fem::AssembledOperator< DF, RF > &op ) {}
          template <class DF, class RF>
          static void end ( Dune::Fem::AssembledOperator< DF, RF > &op ) { op.flushAssembly(); }
        };

      } // namespace Global




      namespace Impl
      {

        template< class LocalMatrix >
        struct LocalMatrixGetter
        {
          typedef typename LocalMatrix::value_type value_type;

          explicit LocalMatrixGetter ( const LocalMatrix &localMatrix ) : localMatrix_( localMatrix ) {}

          decltype( auto ) operator[] ( int row ) const { return localMatrix_[ row ]; }

          value_type get ( int row, int col ) const { return localMatrix_[ row ][ col ]; }

          const LocalMatrix & localMatrix() const { return localMatrix_; }

        private:
          const LocalMatrix &localMatrix_;
        };

      } // namespace Impl



      // AddBase
      // -------

      template< class AssembledOperator >
      struct AddBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
      {
        typedef typename AssembledOperator::RangeFieldType value_type;

        typedef Global::Add< AssembledOperator > GlobalOperationType;

        template< class DomainEntity, class RangeEntity, class LocalMatrix >
        void begin ( const DomainEntity &domainEntity, const RangeEntity &rangeEntity, const AssembledOperator &op, LocalMatrix &localMatrix ) const
        {
          localMatrix.clear();
        }

        template< class DomainEntity, class RangeEntity, class LocalMatrix >
        void end ( const DomainEntity &domainEntity, const RangeEntity &rangeEntity, LocalMatrix &localMatrix, AssembledOperator &op ) const
        {
          op.addLocalMatrix( domainEntity, rangeEntity, Impl::LocalMatrixGetter< LocalMatrix >( localMatrix ) );
        }
      };



      // AddScaledBase
      // -------------

      template< class AssembledOperator >
      struct AddScaledBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
        : public AddBase< AssembledOperator >
      {
        typedef typename AssembledOperator::RangeFieldType value_type;

        AddScaledBase ( value_type factor ) : factor_( std::move( factor ) ) {}

        template< class DomainEntity, class RangeEntity, class LocalMatrix >
        void end ( const DomainEntity &domainEntity, const RangeEntity &rangeEntity, LocalMatrix &localMatrix, AssembledOperator &op ) const
        {
          op.addScaledLocalMatrix( domainEntity, rangeEntity, Impl::LocalMatrixGetter< LocalMatrix >( localMatrix ), factor_ );
        }

      private:
        value_type factor_;
      };



      // SetAndSelectBaseImpl
      // ---------------------

      template< class AssembledOperator, const bool getAndSet >
      struct SetAndSelectBaseImpl
      {
        typedef typename AssembledOperator::RangeFieldType value_type;

        typedef Global::Set< AssembledOperator > GlobalOperationType;

        template< class DomainEntity, class RangeEntity, class LocalMatrix >
        void begin ( const DomainEntity &domainEntity, const RangeEntity &rangeEntity, const AssembledOperator &op, LocalMatrix &localMatrix ) const
        {
          if constexpr ( getAndSet )
          {
            // get values first for modification
            op.getLocalMatrix( domainEntity, rangeEntity, localMatrix );
          }
          else // set only
          {
            localMatrix.clear();
          }
        }

        template< class DomainEntity, class RangeEntity, class LocalMatrix >
        void end ( const DomainEntity &domainEntity, const RangeEntity &rangeEntity, LocalMatrix &localMatrix, AssembledOperator &op ) const
        {
          op.setLocalMatrix( domainEntity, rangeEntity, Impl::LocalMatrixGetter< LocalMatrix >( localMatrix ) );
        }
      };

      // SetBase
      // -------
      template< class AssembledOperator >
      struct SetBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
        : public SetAndSelectBaseImpl< AssembledOperator, false > // set only
      {
      };

      // SetSelectedBase
      // ---------------

      template< class AssembledOperator >
      struct SetSelectedBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
        : public SetAndSelectBaseImpl< AssembledOperator, true > // get and set
      {
      };

    } // namespace Assembly

  } // namespace Fem


  namespace Fem
  {

    // LocalContribution for Assembled Operators
    // -----------------------------------------

    template< class AssembledOperator, template< class > class AssemblyOperation >
    class LocalContribution< AssembledOperator, AssemblyOperation, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
      : public TemporaryLocalMatrix< typename AssembledOperator::DomainFunctionType::DiscreteFunctionSpaceType,
                                     typename AssembledOperator::RangeFunctionType::DiscreteFunctionSpaceType
                                   >
    {
    public:
      typedef AssembledOperator AssembledOperatorType;
      typedef AssemblyOperation< AssembledOperator > AssemblyOperationType;

      typedef typename AssembledOperatorType::DomainFunctionType::DiscreteFunctionSpaceType DomainSpaceType;
      typedef typename AssembledOperatorType::RangeFunctionType::DiscreteFunctionSpaceType  RangeSpaceType;

    private:
      typedef LocalContribution< AssembledOperator, AssemblyOperation > ThisType;
      typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType >   BaseType;

    public:
      typedef typename DomainSpaceType::BasisFunctionSetType   DomainBasisFunctionSetType;
      typedef typename RangeSpaceType::BasisFunctionSetType    RangeBasisFunctionSetType;

      typedef typename DomainBasisFunctionSetType::EntityType  DomainEntityType;
      typedef typename RangeBasisFunctionSetType::EntityType   RangeEntityType;

      typedef typename AssembledOperatorType::RangeFieldType value_type;

      typedef typename BaseType :: MatrixEntriesType  LocalMatrixEntriesType;
      typedef typename LocalMatrixEntriesType::size_type SizeType;

      typedef RowReferenceVector< value_type > row_reference;
      typedef RowReferenceVector< const value_type > const_row_reference;

    private:
      template< class T, class... U >
      struct InTypeRange
        : public std::integral_constant< bool, Std::Or( std::is_same< T, U >::value... ) >
      {};

      template< class T >
      struct IsRangeValue
        : public InTypeRange< T, typename RangeBasisFunctionSetType::RangeType,
                                 typename RangeBasisFunctionSetType::JacobianRangeType,
                                 typename RangeBasisFunctionSetType::HessianRangeType >
      {};

      struct ColIndexMapper
      {
        ColIndexMapper ( SizeType j, SizeType cols ) : j_( j ), cols_( cols ) {}

        SizeType operator[] ( SizeType i ) const { return (i*cols_ + j_); }

      private:
        SizeType j_, cols_;
      };

    protected:
      using BaseType::mat_cols;

    public:
      using BaseType::domainBasisFunctionSet;
      using BaseType::rangeBasisFunctionSet;

      template< class... Args >
      explicit LocalContribution ( AssembledOperator &assembledOperator, Args &&... args )
        : BaseType( assembledOperator.domainSpace(), assembledOperator.rangeSpace() ),
          assembledOperator_( assembledOperator ),
          assemblyOperation_( std::forward< Args >( args )... )
      {
        assembledOperator.template beginAssemble< typename AssemblyOperationType::GlobalOperationType >();
      }

      LocalContribution ( const ThisType & ) = delete;
      LocalContribution ( ThisType && ) = delete;

      ~LocalContribution () { assembledOperator().template endAssemble< typename AssemblyOperationType::GlobalOperationType >(); }

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      const AssembledOperatorType &assembledOperator () const { return assembledOperator_; }
      AssembledOperatorType &assembledOperator () { return assembledOperator_; }

      SubVector< LocalMatrixEntriesType, ColIndexMapper > column ( SizeType j )
      {
        return SubVector< LocalMatrixEntriesType, ColIndexMapper >( localMatrixEntries(), ColIndexMapper( j, mat_cols() ) );
      }

      SubVector< const LocalMatrixEntriesType, ColIndexMapper > column ( SizeType j ) const
      {
        return SubVector< const LocalMatrixEntriesType, ColIndexMapper >( localMatrixEntries(), ColIndexMapper( j, mat_cols() ) );
      }

      // this method behaves different to column
      typename BaseType::MatrixColumnType matrixColumn( SizeType j )
      {
        return BaseType::column( j );
      }

      // inherited from DenseMatrix
      using BaseType :: axpy;

      template< class Point, class... Factors >
      auto axpy ( const Point &x, const Factors &... factors )
        -> std::enable_if_t< Std::And( (IsRangeValue< std::decay_t< decltype( std::declval< Factors & >()[ 0 ] ) > >::value)... ) >
      {
        const SizeType matCols = mat_cols();
        for( SizeType j = 0; j < matCols; ++j )
        {
          auto col = column( j );
          rangeBasisFunctionSet().axpy( x, factors[ j ]..., col );
        }
      }

      /**
       * \brief obtain the order of this local contribution
       *
       * The order of a local contribution refers to the polynomial order required
       * to integrate it exactly.
       *
       * \note It is not completely clear what this value should be, e.g., for
       *       bilinear basis functions.
       *
       * \returns order of the local contribution
       **/
      int order () const { return domainBasisFunctionSet().order() + rangeBasisFunctionSet().order(); }

      /** \copydoc Dune::Fem::LocalMatrixInterface::bind */
      void bind ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        BaseType::bind( domainEntity, rangeEntity );
        assemblyOperation_.begin( domainEntity, rangeEntity, assembledOperator(), *this );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::unbind */
      void unbind ()
      {
        assemblyOperation_.end( domainBasisFunctionSet().entity(), rangeBasisFunctionSet().entity(), *this, assembledOperator() );
        BaseType::unbind();
      }

      /** \brief return const reference to vector of local matrix entries **/
      const LocalMatrixEntriesType &localMatrixEntries () const { return fields_; }
      /** \brief return reference to vector of local matrix entries **/
      LocalMatrixEntriesType &localMatrixEntries () { return fields_; }

    protected:
      using BaseType::fields_;
      AssembledOperatorType &assembledOperator_;
      AssemblyOperationType assemblyOperation_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_LOCALCONTRIBUTION_HH
