#ifndef DUNE_FEM_DENSEMATRIX_HH
#define DUNE_FEM_DENSEMATRIX_HH

#include <dune/fem/storage/objectstack.hh>
#include <dune/fem/solver/oemsolver.hh>
#include <dune/fem/function/adaptivefunction.hh>
#include <dune/fem/misc/functor.hh>
#include <dune/fem/operator/common/localmatrix.hh>
#include <dune/fem/operator/common/localmatrixwrapper.hh>

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

      DenseRowMatrix ()
      : rows_( 0 ),
        cols_( 0 ),
        fields_( 0 )
      {}

      DenseRowMatrix ( unsigned int rows, unsigned int cols )
      : rows_( 0 ),
        cols_( 0 ),
        fields_( 0 )
      {
        reserve( rows, cols );
      }

      ~DenseRowMatrix ()
      {
        delete[] fields_;
      }

      unsigned int rows () const
      {
        return rows_;
      }

      unsigned int cols () const
      {
        return cols_;
      }

      const Field &operator() ( const unsigned int row, const unsigned int col ) const
      {
        assert( (row < rows()) && (col < cols()) );
        return fields_[ row*cols() + col ];
      }

      Field &operator() ( const unsigned int row, const unsigned int col )
      {
        assert( (row < rows()) && (col < cols()) );
        return fields_[ row*cols() + col ];
      }

      void add(const unsigned int row, const unsigned int col, const Field& value )
      {
        assert( (row < rows()) && (col < cols()) );
        fields_[ row*cols() + col ] += value;
      }

      Row< const Field > operator[] ( const unsigned int row ) const
      {
        assert( row < rows() );
        return Row< const Field >( cols(), fields_ + row*cols() );
      }

      Row< Field > operator[] ( const unsigned int row )
      {
        assert( row < rows() );
        return Row< Field >( cols(), fields_ + row*cols() );
      }

      void clear ()
      {
        Field *end = fields_ + (rows_*cols_);
        for( Field *it = fields_; it != end; ++it )
          *it = Field( 0 );
      }

      void multiply ( const Field *x, Field *y ) const
      {
        for( unsigned int row = 0; row < rows(); ++row )
        {
          const Field *fields = fields_ + row*cols();
          y[ row ] = Field( 0 );
          for( unsigned int col = 0; col < cols(); ++col )
            y[ row ] += fields[ col ] * x[ col ];
        }
      }

      void reserve ( unsigned int rows, unsigned int cols )
      {
        if( (rows != rows_) || (cols != cols_) )
        {
          delete[] fields_;
          fields_ = new Field[ rows*cols ];
          rows_ = rows;
          cols_ = cols;
        }
        // Martin: Is clearing required, here?
        clear();
      }

      void print( std::ostream& out ) const 
      {
        for( unsigned int row = 0; row < rows(); ++row )
        {
          const Field *fields = fields_ + row*cols();
          for( unsigned int col = 0; col < cols(); ++col )
            out << fields[ col ] << " ";

          out << std::endl;
        }
      }

    private:
      unsigned int rows_;
      unsigned int cols_;
      Field *fields_;
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
    : public OEMMatrix
    {
      typedef DenseRowMatrixObject< DomainSpace, RangeSpace > ThisType;

    public:
      typedef DomainSpace DomainSpaceType;
      typedef RangeSpace RangeSpaceType;

      typedef typename RangeSpaceType::RangeFieldType Field;

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

      const DomainSpace &domainSpace () const
      {
        return domainSpace_;
      }

      const RangeSpace &rangeSpace () const
      {
        return rangeSpace_;
      }

    private:
      const DomainSpaceType &domainSpace_;
      const RangeSpaceType &rangeSpace_;

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

      LocalMatrix ( MatrixType &matrix,
                    const DomainSpaceType &domainSpace,
                    const RangeSpaceType &rangeSpace )
      : BaseType( domainSpace, rangeSpace ),
        matrix_( matrix )
      {}

    private:
      LocalMatrix ( const ThisType & );
      ThisType &operator= ( const ThisType & );

    public:
      void init ( const RowEntityType &rowEntity, const ColEntityType &colEntity )
      {
        BaseType::init( rowEntity, colEntity );
        
        map( rangeSpace().mapper(), rowEntity, rowIndices_ );
        map( domainSpace().mapper(), colEntity, colIndices_ );
      }

      int rows () const
      {
        return rowIndices_.size();
      }

      int cols () const
      {
        return colIndices_.size();
      }

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

    private:
      template< class Mapper, class Entity >
      void map ( const Mapper &mapper, const Entity &entity, std::vector< unsigned int > &indices )
      {
        indices.resize( mapper.numDofs( entity ) );
        mapper.mapEach( entity, Fem::AssignFunctor< std::vector< unsigned int > >( indices ) );
      }

    protected:
      using BaseType::domainSpace_;
      using BaseType::rangeSpace_;

    private:
      MatrixType &matrix_;

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
        return new ObjectType( matrixObject_->matrix_, matrixObject_->domainSpace_, matrixObject_->rangeSpace_ );
      }

    private:
      MatrixObject *matrixObject_;
    };

  } // namespace Fem

} // namespace Dune


// the following is deprecated
#include <dune/fem/operator/common/stencil.hh>

namespace Dune
{

  namespace Fem
  {

    // DenseRowMatrixOperator
    // ----------------------

    template< class DomainFunction, class RangeFunction , class Traits >
    class DenseRowMatrixOperator
    : public DenseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType >,
      public AssembledOperator< DomainFunction, RangeFunction >
    {
      typedef DenseRowMatrixOperator< DomainFunction, RangeFunction, Traits > ThisType;
      typedef DenseRowMatrixObject< typename DomainFunction::DiscreteFunctionSpaceType, typename RangeFunction::DiscreteFunctionSpaceType > BaseType;

    public:
      typedef typename BaseType::DomainSpaceType DomainSpaceType;
      typedef typename BaseType::RangeSpaceType RangeSpaceType;

      using BaseType::apply;

      DUNE_VERSION_DEPRECATED(1,4,remove)
      DenseRowMatrixOperator ( const std::string &name,
                               const DomainSpaceType &domainSpace,
                               const RangeSpaceType &rangeSpace )
      : BaseType( domainSpace, rangeSpace )
      {}

      virtual void operator() ( const DomainFunction &u, RangeFunction &w ) const
      {
        apply( u, w );
      }

      const BaseType &systemMatrix () const { return *this; }

      void communicate () const {}

      void reserve ( bool verbose = false )
      {
        BaseType::reserve( DummyStencil(), verbose );
      }

    private:
      class DummyStencil {};
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DENSEMATRIX_HH
