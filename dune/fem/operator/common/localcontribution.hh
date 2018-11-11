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
          static void begin ( AssembledOperator &op ) {}
          static void end ( AssembledOperator &op ) { op.communicate(); }
        };



        // SetBase
        // -------

        template< class AssembledOperator >
        struct SetBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
        {
          static void begin ( AssembledOperator &op ) {}
          static void end ( AssembledOperator &op ) { op.communicate(); }
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



      // SetBase
      // -------

      template< class AssembledOperator >
      struct SetBase< AssembledOperator, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
      {
        typedef typename AssembledOperator::RangeFieldType value_type;

        typedef Global::Set< AssembledOperator > GlobalOperationType;

        template< class DomainEntity, class RangeEntity, class LocalMatrix >
        void begin ( const DomainEntity &domainEntity, const RangeEntity &rangeEntity, const AssembledOperator &op, LocalMatrix &localMatrix ) const
        {
          localMatrix.clear();
        }

        template< class DomainEntity, class RangeEntity, class LocalMatrix >
        void end ( const DomainEntity &domainEntity, const RangeEntity &rangeEntity, LocalMatrix &localMatrix, AssembledOperator &op ) const
        {
          op.setLocalMatrix( domainEntity, rangeEntity, Impl::LocalMatrixGetter< LocalMatrix >( localMatrix ) );
        }
      };

    } // namespace Assembly

  } // namespace Fem



  // DenseMatVecTraits for LocalContribution
  // ---------------------------------------

  template< class AssembledOperator, template< class > class AssemblyOperation >
  struct DenseMatVecTraits< Fem::LocalContribution< AssembledOperator, AssemblyOperation, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > > >
  {
    typedef Fem::LocalContribution< AssembledOperator, AssemblyOperation > derived_type;
    typedef typename AssembledOperator::RangeFieldType value_type;
    typedef DynamicVector< value_type > row_type;
    typedef std::vector< value_type > container_type;
    typedef typename container_type::size_type size_type;

    typedef Fem::RowReferenceVector< value_type > row_reference;
    typedef Fem::RowReferenceVector< const value_type > const_row_reference;
  };



  // FieldTraits for LocalContribution
  // ---------------------------------

  template< class AssembledOperator, template< class > class AssemblyOperation >
  struct FieldTraits< Fem::LocalContribution< AssembledOperator, AssemblyOperation, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > > >
  {
    typedef typename FieldTraits< typename AssembledOperator::RangeFieldType >::field_type field_type;
    typedef typename FieldTraits< typename AssembledOperator::RangeFieldType >::real_type real_type;
  };



  namespace Fem
  {

    // LocalContribution for Assembled Operators
    // -----------------------------------------

    template< class AssembledOperator, template< class > class AssemblyOperation >
    class LocalContribution< AssembledOperator, AssemblyOperation, std::enable_if_t< Fem::IsAssembledOperator< AssembledOperator >::value > >
      : public Dune::DenseMatrix< LocalContribution< AssembledOperator, AssemblyOperation > >
    {
      typedef LocalContribution< AssembledOperator, AssemblyOperation > ThisType;
      typedef Dune::DenseMatrix< LocalContribution< AssembledOperator, AssemblyOperation > > BaseType;

    public:
      typedef AssembledOperator AssembledOperatorType;
      typedef AssemblyOperation< AssembledOperator > AssemblyOperationType;

      typedef typename AssembledOperatorType::DomainFunctionType::DiscreteFunctionSpaceType::BasisFunctionSetType DomainBasisFunctionSetType;
      typedef typename AssembledOperatorType::RangeFunctionType::DiscreteFunctionSpaceType::BasisFunctionSetType RangeBasisFunctionSetType;

      typedef typename DomainBasisFunctionSetType::EntityType DomainEntityType;
      typedef typename RangeBasisFunctionSetType::EntityType RangeEntityType;

      typedef typename AssembledOperatorType::RangeFieldType value_type;

      typedef std::vector< value_type > LocalMatrixEntriesType;
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
        : public InTypeRange< T, typename RangeBasisFunctionSetType::RangeType, typename RangeBasisFunctionSetType::JacobianRangeType, typename RangeBasisFunctionSetType::HessianRangeType >
      {};

      struct ColIndexMapper
      {
        ColIndexMapper ( SizeType j, SizeType cols ) : j_( j ), cols_( cols ) {}

        SizeType operator[] ( SizeType i ) const { return (i*cols_ + j_); }

      private:
        SizeType j_, cols_;
      };

    public:
      template< class... Args >
      explicit LocalContribution ( AssembledOperator &assembledOperator, Args &&... args )
        : assembledOperator_( assembledOperator ),
          localMatrixEntries_( assembledOperator.domainSpace().maxNumDofs() * assembledOperator.rangeSpace().maxNumDofs() ),
          assemblyOperation_( std::forward< Args >( args )... )
      {
        //assembledOperator.template beginAssemble< typename AssemblyOperationType::GlobalOperationType >();
      }

      LocalContribution ( const ThisType & ) = delete;
      LocalContribution ( ThisType && ) = delete;

      ~LocalContribution () { /*assembledOperator().template endAssemble< typename AssemblyOperationType::GlobalOperationType >();*/ }

      ThisType &operator= ( const ThisType & ) = delete;
      ThisType &operator= ( ThisType && ) = delete;

      const AssembledOperatorType &assembledOperator () const { return assembledOperator_; }
      AssembledOperatorType &assembledOperator () { return assembledOperator_; }

      /** \brief set all DoFs to zero **/
      void clear () { std::fill( localMatrixEntries().begin(), localMatrixEntries().end(), value_type( 0 ) ); }

      SubVector< LocalMatrixEntriesType, ColIndexMapper > column ( SizeType j )
      {
        return SubVector< LocalMatrixEntriesType, ColIndexMapper >( localMatrixEntries(), ColIndexMapper( j, mat_cols() ) );
      }

      SubVector< const LocalMatrixEntriesType, ColIndexMapper > column ( SizeType j ) const
      {
        return SubVector< const LocalMatrixEntriesType, ColIndexMapper >( localMatrixEntries(), ColIndexMapper( j, mat_cols() ) );
      }

      template< class Point, class... Factors >
      auto axpy ( const Point &x, const Factors &... factors )
        -> std::enable_if_t< Std::And( (IsRangeValue< std::decay_t< decltype( std::declval< Factors & >()[ 0 ] ) > >::value)... ) >
      {
        for( SizeType j = 0; j < mat_cols(); ++j )
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

      /**
       * \brief obtain the basis function set for the domain of this local contribution
       * \returns reference to the basis function set
       **/
      const DomainBasisFunctionSetType &domainBasisFunctionSet () const { return domainBasisFunctionSet_; }

      /**
       * \brief obtain the basis function set for the range of this local contribution
       * \returns reference to the basis function set
       **/
      const RangeBasisFunctionSetType &rangeBasisFunctionSet () const { return rangeBasisFunctionSet_; }

      void bind ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        domainBasisFunctionSet_ = assembledOperator().domainSpace().basisFunctionSet( domainEntity );
        rangeBasisFunctionSet_ = assembledOperator().rangeSpace().basisFunctionSet( rangeEntity );
        localMatrixEntries().resize( mat_rows() * mat_cols() );
        assemblyOperation_.begin( domainEntity, rangeEntity, assembledOperator(), *this );
      }

      void unbind ()
      {
        assemblyOperation_.end( domainBasisFunctionSet().entity(), rangeBasisFunctionSet().entity(), *this, assembledOperator() );
        domainBasisFunctionSet_ = DomainBasisFunctionSetType();
        rangeBasisFunctionSet_ = RangeBasisFunctionSetType();
      }

      /** \brief return const reference to vector of local matrix entries **/
      const LocalMatrixEntriesType &localMatrixEntries () const { return localMatrixEntries_; }
      /** \brief return reference to vector of local matrix entries **/
      LocalMatrixEntriesType &localMatrixEntries () { return localMatrixEntries_; }

      // make this thing a dense matrix
      SizeType mat_rows () const { return rangeBasisFunctionSet().size(); }
      SizeType mat_cols () const { return domainBasisFunctionSet().size(); }

      row_reference mat_access ( SizeType i )
      {
        const SizeType cols = mat_cols();
        return row_reference( localMatrixEntries_.data() + i*cols, cols );
      }

      const_row_reference mat_access ( SizeType i ) const
      {
        const SizeType cols = mat_cols();
        return const_row_reference( localMatrixEntries_.data() + i*cols, cols );
      }

      const value_type* data() const { return localMatrixEntries_.data(); }

    private:
      AssembledOperatorType &assembledOperator_;
      LocalMatrixEntriesType localMatrixEntries_;
      DomainBasisFunctionSetType domainBasisFunctionSet_;
      RangeBasisFunctionSetType rangeBasisFunctionSet_;
      AssemblyOperationType assemblyOperation_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_COMMON_LOCALCONTRIBUTION_HH
