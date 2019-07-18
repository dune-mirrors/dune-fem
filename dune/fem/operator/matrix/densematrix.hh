#ifndef DUNE_FEM_OPERATOR_MATRIX_DENSEMATRIX_HH
#define DUNE_FEM_OPERATOR_MATRIX_DENSEMATRIX_HH

#include <memory>
#include <iostream>

#include <dune/common/dynmatrixev.hh>

#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>
#include <dune/fem/solver/krylovinverseoperators.hh>
#include <dune/fem/storage/objectstack.hh>

namespace Dune
{

  namespace Fem
  {

    // DenseRowMatrix
    // --------------

    template< class F >
    class DenseRowMatrix
    {
      typedef DenseRowMatrix< F > ThisType;

    public:
      typedef F Field;

      template< class RF >
      class Row;

      DenseRowMatrix () : rows_( 0 ), cols_( 0 ) {}

      DenseRowMatrix ( unsigned int rows, unsigned int cols )
        : rows_( 0 ), cols_( 0 )
      {
        reserve( rows, cols );
      }

      unsigned int rows () const { return rows_; }
      unsigned int cols () const { return cols_; }

      const Field &operator() ( unsigned int row, unsigned int col ) const
      {
        assert( (row < rows()) && (col < cols()) );
        return fields_[ row*cols() + col ];
      }

      Field &operator() ( unsigned int row, unsigned int col )
      {
        assert( (row < rows()) && (col < cols()) );
        return fields_[ row*cols() + col ];
      }

      void add ( unsigned int row, unsigned int col, const Field &value )
      {
        assert( (row < rows()) && (col < cols()) );
        fields_[ row*cols() + col ] += value;
      }

      Row< const Field > operator[] ( unsigned int row ) const
      {
        assert( row < rows() );
        return Row< const Field >( cols(), fields_.get() + row*cols() );
      }

      Row< Field > operator[] ( unsigned int row )
      {
        assert( row < rows() );
        return Row< Field >( cols(), fields_.get() + row*cols() );
      }

      void clear () { std::fill( fields_.get(), fields_.get() + (rows()*cols()), Field( 0 ) ); }

      void multiply ( const Field *x, Field *y ) const
      {
        for( unsigned int row = 0; row < rows(); ++row )
        {
          const Field *fields = fields_.get() + row*cols();
          y[ row ] = Field( 0 );
          for( unsigned int col = 0; col < cols(); ++col )
            y[ row ] += fields[ col ] * x[ col ];
        }
      }

      /**
       * \brief calculate eigenvalues
       *
       * \returns  eigen values in ascending order
       *
       * \note LAPACK::dgeev is used to compute the eigen values.
       *
       * \note The matrix is destroyed.
       */
      std::vector< std::complex< double > > eigenValues ()
      {
        if( rows() != cols() )
          DUNE_THROW( InvalidStateException, "Requiring square matrix for eigenvalue computation" );

        const long int N = rows();
        const char jobvl = 'n';
        const char jobvr = 'n';

        // working memory
        std::unique_ptr< double[] > work = std::make_unique< double[] >( 5*N );

        // return value information
        long int info = 0;
        long int lwork = 3*N;

        // call LAPACK routine (see fmatrixev_ext.cc)
        DynamicMatrixHelp::eigenValuesNonsymLapackCall( &jobvl, &jobvr, &N, fields_, &N, work.get(), work.get()+N, 0, &N, 0, &N, work.get()+2*N, &lwork, &info );

        if( info != 0 )
          DUNE_THROW( MathError, "DenseRowMatrix: Eigenvalue computation failed" );

        std::vector< std::complex< double > > eigenValues( N );
        std::transform( work.get(), work.get()+N, work.get()+N, eigenValues.begin(), [] ( double r, double i ) { return std::complex< double >( r, i ); } );
        return eigenValues;
      }

      void reserve ( unsigned int rows, unsigned int cols )
      {
        if( (rows != rows_) || (cols != cols_) )
        {
          fields_.reset( new Field[ rows*cols ] );
          rows_ = rows;
          cols_ = cols;
        }
        // Martin: Is clearing required, here?
        clear();
      }

      void print( std::ostream& s=std::cout ) const
      {
        s.precision( 6 );
        for( unsigned int row = 0; row < rows(); ++row )
        {
          const Field *fields = fields_ + row*cols();
          for( unsigned int col = 0; col < cols(); ++col )
            s << fields[ col ] << " ";

          s << std::endl;
        }
      }

    private:
      unsigned int rows_, cols_;
      std::unique_ptr< Field[] > fields_;
    };



    // DenseRowMatrix::Row
    // -------------------

    template< class F >
    template< class RF >
    class DenseRowMatrix< F >::Row
    {
      typedef Row< RF > ThisType;

      template< class > friend class Row;

    public:
      Row ( unsigned int cols, RF *fields )
      : cols_( cols ),
        fields_( fields )
      {}

      Row ( const Row< F > &row )
      : cols_( row.cols_ ),
        fields_( row.fields_ )
      {}

      const RF &operator[] ( const unsigned int col ) const
      {
        assert( col < size() );
        return fields_[ col ];
      }

      RF &operator[] ( const unsigned int col )
      {
        assert( col < size() );
        return fields_[ col ];
      }

      void clear ()
      {
        Field *const end = fields_ + size();
        for( Field *it = fields_; it != end; ++it )
          *it = Field( 0 );
      }

      unsigned int size () const
      {
        return cols_;
      }

    private:
      unsigned int cols_;
      RF *fields_;
    };



    // DenseRowMatrixObject
    // --------------------

    template< class DomainSpace, class RangeSpace >
    class DenseRowMatrixObject
    {
      typedef DenseRowMatrixObject< DomainSpace, RangeSpace > ThisType;

    public:
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;

      typedef typename RangeSpaceType::RangeFieldType Field;

      typedef typename DomainSpaceType::BlockMapperType DomainBlockMapperType;
      typedef NonBlockMapper< DomainBlockMapperType, DomainSpaceType::localBlockSize > DomainMapperType;
      typedef typename RangeSpaceType::BlockMapperType RangeBlockMapperType;
      typedef NonBlockMapper< RangeBlockMapperType, RangeSpaceType::localBlockSize > RangeMapperType;

      typedef typename DomainSpaceType::EntityType DomainEntityType;
      typedef typename RangeSpaceType::EntityType RangeEntityType;
      typedef typename DomainSpace::GridType::template Codim< 0 >::Entity ColEntityType;
      typedef typename RangeSpace::GridType::template Codim< 0 >::Entity RowEntityType;

      typedef DenseRowMatrix< Field > MatrixType;

    private:
      class LocalMatrixTraits;
      class LocalMatrix;
      class LocalMatrixFactory;

      typedef Fem::ObjectStack< LocalMatrixFactory > LocalMatrixStack;

    public:
      typedef LocalMatrixWrapper< LocalMatrixStack > LocalMatrixType;

      DenseRowMatrixObject ( const DomainSpaceType &domainSpace,
                             const RangeSpaceType &rangeSpace )
      : domainSpace_( domainSpace ),
        rangeSpace_( rangeSpace ),
        domainMapper_( domainSpace.blockMapper() ),
        rangeMapper_( rangeSpace.blockMapper() ),
        domainSequence_( -1 ),
        rangeSequence_( -1 ),
        localMatrixFactory_( *this ),
        localMatrixStack_( localMatrixFactory_ )
      {}

      MatrixType &matrix ()
      {
        return matrix_;
      }

      LocalMatrixType localMatrix ( const RowEntityType &rowEntity,
                                    const ColEntityType &colEntity )
      {
        return LocalMatrixType( localMatrixStack_, rowEntity, colEntity );
      }


      LocalMatrixType localMatrix () const
      {
        return LocalMatrixType( localMatrixStack_ );
      }

      void clear ()
      {
        matrix_.clear();
      }

      template< class Stencil >
      void reserve ( const Stencil &stencil, bool verbose = false )
      {
        if( (domainSequence_ != domainSpace().sequence()) || (rangeSequence_ != rangeSpace().sequence()) )
        {
          matrix_.reserve( rangeSpace().size(), domainSpace().size() );
          domainSequence_ = domainSpace().sequence();
          rangeSequence_ = rangeSpace().sequence();
        }
      }

      template< class DomainFunction, class RangeFunction >
      void apply ( const DomainFunction &u, RangeFunction &w ) const
      {
        matrix_.multiply( u.leakPointer(), w.leakPointer() );
        rangeSpace().communicate( w );
      }

      Field ddotOEM ( const Field *v, const Field *w ) const
      {
        typedef AdaptiveDiscreteFunction< RangeSpaceType > RangeFunction;
        RangeFunction vFunction( "v (ddotOEM)", rangeSpace(), v );
        RangeFunction wFunction( "w (ddotOEM)", rangeSpace(), w );
        return vFunction.scalarProductDofs( wFunction );
      }

      void multOEM ( const Field *u, Field *w ) const
      {
        matrix_.multiply( u, w );

        typedef AdaptiveDiscreteFunction< RangeSpaceType > RangeFunction;
        RangeFunction wFunction( "w (multOEM)", rangeSpace(), w );
        rangeSpace().communicate( wFunction );
      }

      template< class DiscreteFunctionType >
      void extractDiagonal( DiscreteFunctionType &diag ) const
      {
        assert( matrix_.rows() == matrix_.cols() );
        const auto dofEnd = diag.dend();
        unsigned int row = 0;
        for( auto dofIt = diag.dbegin(); dofIt != dofEnd; ++dofIt, ++row )
        {
          assert( row < matrix_.rows() );
          (*dofIt) = matrix_( row, row );
        }
      }

      const DomainSpace &domainSpace () const { return domainSpace_; }
      const RangeSpace &rangeSpace () const { return rangeSpace_; }

    private:
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;
      DomainMapperType domainMapper_;
      RangeMapperType rangeMapper_;

      int domainSequence_;
      int rangeSequence_;

      MatrixType matrix_;

      LocalMatrixFactory localMatrixFactory_;
      mutable LocalMatrixStack localMatrixStack_;
    };



    // DenseRowMatrixObject::LocalMatrixTraits
    // ---------------------------------------

    template< class DomainSpace, class RangeSpace >
    class DenseRowMatrixObject< DomainSpace, RangeSpace >::LocalMatrixTraits
    {
      typedef DenseRowMatrixObject< DomainSpace, RangeSpace > MatrixObject;

    public:
      typedef typename MatrixObject::DomainSpaceType DomainSpaceType;
      typedef typename MatrixObject::RangeSpaceType RangeSpaceType;

      typedef typename MatrixObject::Field RangeFieldType;
      typedef RangeFieldType LittleBlockType;

      typedef typename MatrixObject::LocalMatrix LocalMatrixType;
    };



    // DenseRowMatrixObject::LocalMatrix
    // ---------------------------------

    template< class DomainSpace, class RangeSpace >
    class DenseRowMatrixObject< DomainSpace, RangeSpace >::LocalMatrix
    : public LocalMatrixDefault< LocalMatrixTraits >
    {
      typedef DenseRowMatrixObject< DomainSpace, RangeSpace > MatrixObject;

      typedef LocalMatrix ThisType;
      typedef LocalMatrixDefault< LocalMatrixTraits > BaseType;

    public:
      typedef LocalMatrixTraits Traits;

      typedef typename MatrixObject::MatrixType MatrixType;

      typedef typename Traits::RangeFieldType RangeFieldType;
      typedef typename Traits::LittleBlockType LittleBlockType;
      typedef RangeFieldType DofType;

      LocalMatrix ( MatrixType &matrix, const DomainSpaceType &domainSpace, const RangeSpaceType &rangeSpace, const DomainMapperType &domainMapper, const RangeMapperType &rangeMapper )
        : BaseType( domainSpace, rangeSpace ),
          matrix_( matrix ),
          domainMapper_( domainMapper ),
          rangeMapper_( rangeMapper )
      {}

      LocalMatrix ( const ThisType & ) = delete;
      ThisType &operator= ( const ThisType & ) = delete;

      void init ( const DomainEntityType &domainEntity, const RangeEntityType &rangeEntity )
      {
        BaseType::init( domainEntity, rangeEntity );

        rowIndices_.resize( rangeMapper_.numDofs( rangeEntity ) );
        rangeMapper_.mapEach( rangeEntity, Fem::AssignFunctor< std::vector< unsigned int > >( rowIndices_ ) );
        colIndices_.resize( domainMapper_.numDofs( domainEntity ) );
        domainMapper_.mapEach( domainEntity, Fem::AssignFunctor< std::vector< unsigned int > >( colIndices_ ) );
      }

      int rows () const { return rowIndices_.size(); }
      int cols () const { return colIndices_.size(); }

      void add ( const int row, const int col, const DofType &value )
      {
        assert( (row >= 0) && (row < rows()) );
        assert( (col >= 0) && (col < cols()) );
        matrix_( rowIndices_[ row ], colIndices_[ col ] ) += value;
      }

      const DofType &get ( const int row, const int col ) const
      {
        assert( (row >= 0) && (row < rows()) );
        assert( (col >= 0) && (col < cols()) );
        return matrix_( rowIndices_[ row ], colIndices_[ col ] );
      }

      void set ( const int row, const int col, const DofType &value )
      {
        assert( (row >= 0) && (row < rows()) );
        assert( (col >= 0) && (col < cols()) );
        matrix_( rowIndices_[ row ], colIndices_[ col ] ) = value;
      }

      void clearRow ( const int row )
      {
        assert( (row >= 0) && (row < rows()) );
        const unsigned int rowIndex = rowIndices_[ row ];
        matrix_[ rowIndex ].clear();
      }

      void unitRow ( const int row )
      {
        clearRow( row );
        set( row, row, DofType( 1 ) );
      }

      void clear ()
      {
        typedef std::vector< unsigned int >::const_iterator Iterator;
        const Iterator rowEnd = rowIndices_.end();
        for( Iterator rowIt = rowIndices_.begin(); rowIt != rowEnd; ++rowIt )
        {
          const Iterator colEnd = colIndices_.end();
          for( Iterator colIt = colIndices_.begin(); colIt != colEnd; ++colIt )
            matrix_( *rowIt, *colIt ) = DofType( 0 );
        }
      }

      void scale ( const DofType &value )
      {
        typedef std::vector< unsigned int >::const_iterator Iterator;
        const Iterator rowEnd = rowIndices_.end();
        for( Iterator rowIt = rowIndices_.begin(); rowIt != rowEnd; ++rowIt )
        {
          const Iterator colEnd = colIndices_.end();
          for( Iterator colIt = colIndices_.begin(); colIt != colEnd; ++colIt )
            matrix_( *rowIt, *colIt ) *= value;
        }
      }

    protected:
      using BaseType::domainSpace_;
      using BaseType::rangeSpace_;

    private:
      MatrixType &matrix_;
      const DomainMapperType &domainMapper_;
      const RangeMapperType &rangeMapper_;
      std::vector< unsigned int > rowIndices_;
      std::vector< unsigned int > colIndices_;
    };



    // DenseRowMatrixObject::LocalMatrixFactory
    // ----------------------------------------

    template< class DomainSpace, class RangeSpace >
    class DenseRowMatrixObject< DomainSpace, RangeSpace >::LocalMatrixFactory
    {
      typedef DenseRowMatrixObject< DomainSpace, RangeSpace > MatrixObject;

    public:
      typedef LocalMatrix ObjectType;

      LocalMatrixFactory ( MatrixObject &matrixObject )
        : matrixObject_( &matrixObject )
      {}

      ObjectType *newObject () const
      {
        return new ObjectType( matrixObject_->matrix_, matrixObject_->domainSpace_, matrixObject_->rangeSpace_, matrixObject_->domainMapper_, matrixObject_->rangeMapper_ );
      }

    private:
      MatrixObject *matrixObject_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_OPERATOR_MATRIX_DENSEMATRIX_HH
