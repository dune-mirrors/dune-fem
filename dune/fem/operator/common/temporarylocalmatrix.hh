#ifndef DUNE_FEM_TEMPORARYLOCALMATRIX_HH
#define DUNE_FEM_TEMPORARYLOCALMATRIX_HH

#include <algorithm>
#include <vector>

#include <dune/common/densematrix.hh>
#include <dune/common/dynvector.hh>

#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/storage/rowreferencevector.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal Forward Declarations
    // -----------------------------

    template< class DomainSpaceImp, class RangeSpaceImp >
    class TemporaryLocalMatrix;



    // TemporaryLocalMatrixTraits
    // --------------------------

    template< class DomainSpaceImp, class RangeSpaceImp >
    struct TemporaryLocalMatrixTraits
    {
      typedef DomainSpaceImp DomainSpaceType;
      typedef RangeSpaceImp RangeSpaceType;

      typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType >
        LocalMatrixType;

      typedef typename DomainSpaceType :: RangeFieldType DomainFieldType;
      typedef typename RangeSpaceType :: RangeFieldType RangeFieldType;
      typedef RangeFieldType LittleBlockType;
    };

  } // namespace Fem



  // DenseMatVecTraits for TemporaryLocalMatrix
  // ------------------------------------------

  template< class DomainSpaceImp, class RangeSpaceImp >
  struct DenseMatVecTraits< Fem::TemporaryLocalMatrix< DomainSpaceImp, RangeSpaceImp > >
  {
    typedef Fem::TemporaryLocalMatrix< DomainSpaceImp, RangeSpaceImp > derived_type;

    typedef typename Fem::TemporaryLocalMatrixTraits< DomainSpaceImp, RangeSpaceImp >::RangeFieldType value_type;
    typedef int size_type;

    typedef DynamicVector< value_type > row_type;

    typedef Fem::RowReferenceVector< value_type > row_reference;
    typedef Fem::RowReferenceVector< const value_type > const_row_reference;
  };



  // FieldTraits for TemporaryLocalMatrix
  // ------------------------------------

  template< class DomainSpaceImp, class RangeSpaceImp >
  struct FieldTraits< Fem::TemporaryLocalMatrix< DomainSpaceImp, RangeSpaceImp > >
    : public FieldTraits< typename Fem::TemporaryLocalMatrixTraits< DomainSpaceImp, RangeSpaceImp >::RangeFieldType >
  {};



  namespace Fem
  {

    // TemporaryLocalMatrix
    // --------------------

    /** \ingroup Matrix
     *  \class TemporaryLocalMatrix
     *  \brief A local matrix with a small array as storage
     *
     *  A TemporaryLocalMatrix is an implementation of the LocalMatrixInterface
     *  storing the matrix values in an array. It is useful when generating
     *  multiple local matrices that shall then be added together.
     *
     *  \note Due to the backing array, accesses to the matrix should be very fast.
     *
     *  \param DomainSpaceImp  DiscreteFunctionSpace modelling the domain
     *  \param RangeSpaceImp   DiscreteFunctionSpace modelling the range
     */
    template< class DomainSpaceImp, class RangeSpaceImp >
    class TemporaryLocalMatrix
      : public DenseMatrix< TemporaryLocalMatrix< DomainSpaceImp, RangeSpaceImp > >,
        public LocalMatrixDefault< TemporaryLocalMatrixTraits< DomainSpaceImp, RangeSpaceImp > >
    {
    public:
      typedef DomainSpaceImp DomainSpaceType;
      typedef RangeSpaceImp RangeSpaceType;

      typedef TemporaryLocalMatrixTraits< DomainSpaceType, RangeSpaceType >
        Traits;

    private:
      typedef TemporaryLocalMatrix< DomainSpaceType, RangeSpaceType >
        ThisType;
      typedef LocalMatrixDefault< Traits > BaseType;

    public:
      typedef typename Traits :: DomainFieldType DomainFieldType;
      typedef typename Traits :: RangeFieldType RangeFieldType;

      typedef int size_type;
      typedef RangeFieldType value_type;

      typedef RowReferenceVector< value_type > row_reference;
      typedef RowReferenceVector< const value_type > const_row_reference;

      typedef std::vector< RangeFieldType > MatrixEntriesType;
    protected:
      MatrixEntriesType fields_;

    public:
      using BaseType::domainBasisFunctionSet;
      using BaseType::rangeBasisFunctionSet;

      inline TemporaryLocalMatrix ( const DomainSpaceType &domainSpace,
                                    const RangeSpaceType &rangeSpace )
      : BaseType( domainSpace, rangeSpace ),
        fields_()
      {}

      template< class DomainEntityType, class RangeEntityType >
      inline TemporaryLocalMatrix ( const DomainSpaceType &domainSpace,
                                    const RangeSpaceType &rangeSpace,
                                    const DomainEntityType &domainEntity,
                                    const RangeEntityType &rangeEntity )
      : BaseType( domainSpace, rangeSpace, domainEntity, rangeEntity ),
        fields_( rows() * columns() )
      {}

      /** \copydoc Dune::Fem::LocalMatrixInterface::init */
      template< class DomainEntityType, class RangeEntityType >
      inline void init ( const DomainEntityType &domainEntity,
                         const RangeEntityType &rangeEntity )
      {
        bind( domainEntity, rangeEntity );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::bind */
      template< class DomainEntityType, class RangeEntityType >
      inline void bind ( const DomainEntityType &domainEntity,
                         const RangeEntityType &rangeEntity )
      {
        BaseType :: bind( domainEntity, rangeEntity );
        fields_.resize( rows() * columns() );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::unbind */
      inline void unbind ()
      {
        BaseType::unbind();
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::add */
      inline void add ( const int localRow,
                        const int localCol,
                        const RangeFieldType &value )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );
        fields_[ localRow * columns() + localCol ] += value;
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::set */
      inline void set ( const int localRow,
                        const int localCol,
                        const RangeFieldType &value )
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );
        fields_[ localRow * columns() + localCol ] = value;
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::get */
      inline const RangeFieldType get ( const int localRow,
                                        const int localCol ) const
      {
        assert( (localRow >= 0) && (localRow < rows()) );
        assert( (localCol >= 0) && (localCol < columns()) );
        return fields_[ localRow * columns() + localCol ];
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::scale */
      inline void scale( const RangeFieldType &value )
      {
        const std::size_t size = fields_.size();
        for( std::size_t i=0; i<size; ++i )
        {
          fields_[ i ] *= value;
        }
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clear */
      inline void clear ()
      {
        std::fill( fields_.begin() , fields_.end() , 0 );
      }

      /** \copydoc Dune::Fem::LocalMatrixInterface::clearRow */
      void clearRow( const int localRow )
      {
        const int col = columns();
        auto start = fields_.begin() + localRow * col;
        auto end   = start + col;
        std::fill( start, end, 0 );
      }

      size_type rows () const { return mat_rows(); }
      size_type cols () const { return mat_cols(); }
      size_type columns () const { return mat_cols(); }

      // make this thing a dense matrix
      size_type mat_rows () const { return rangeBasisFunctionSet().size(); }
      size_type mat_cols () const { return domainBasisFunctionSet().size(); }

      row_reference mat_access ( size_type i )
      {
        const size_type cols = mat_cols();
        return row_reference( fields_.data() + i*cols, cols );
      }

      const_row_reference mat_access ( size_type i ) const
      {
        const size_type cols = mat_cols();
        return const_row_reference( fields_.data() + i*cols, cols );
      }

      // return pointer to data array for PetscLinearOperator to avoid copying.
      const RangeFieldType* data() const { return fields_.data(); }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_TEMPORARYLOCALMATRIX_HH
