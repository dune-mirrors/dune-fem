#ifndef DUNE_FEM_TEMPORARYLOCALMATRIX_HH
#define DUNE_FEM_TEMPORARYLOCALMATRIX_HH

#include <dune/fem/storage/array.hh>
#include <dune/fem/operator/common/localmatrix.hh>

namespace Dune
{

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
  class TemporaryLocalMatrix;


  
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


  
  template< class DomainSpaceImp, class RangeSpaceImp >
  class TemporaryLocalMatrix
  : public LocalMatrixDefault
    < TemporaryLocalMatrixTraits< DomainSpaceImp, RangeSpaceImp > >
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
    using BaseType :: rows;
    using BaseType :: columns;

  public:
    typedef typename Traits :: DomainFieldType DomainFieldType;
    typedef typename Traits :: RangeFieldType RangeFieldType;

  protected:
    DynamicArray< RangeFieldType > fields_;

  public:
    inline TemporaryLocalMatrix ( const DomainSpaceType &domainSpace,
                                  const RangeSpaceType &rangeSpace )
    : BaseType( domainSpace, rangeSpace ),
      fields_()
    {
    }

    template< class DomainEntityType, class RangeEntityType >
    inline TemporaryLocalMatrix ( const DomainSpaceType &domainSpace,
                                  const RangeSpaceType &rangeSpace,
                                  const DomainEntityType &domainEntity,
                                  const RangeEntityType &rangeEntity )
    : BaseType( domainSpace, rangeSpace, domainEntity, rangeEntity ),
      fields_( rows() * columns() )
    {
    }

    
    /** \copydoc Dune::LocalMatrixInterface::init */
    template< class DomainEntityType, class RangeEntityType >
    inline void init ( const DomainEntityType &domainEntity,
                       const RangeEntityType &rangeEntity )
    {
      BaseType :: init( domainEntity, rangeEntity );
      fields_.resize( rows() * columns() );
    }

    /** \copydoc Dune::LocalMatrixInterface::add */
    inline void add ( const int localRow,
                      const int localCol,
                      const RangeFieldType &value )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );
      fields_[ localRow * columns() + localCol ] += value;
    }

    /** \copydoc Dune::LocalMatrixInterface::set */
    inline void set ( const int localRow, 
                      const int localCol, 
                      const RangeFieldType &value )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );
      fields_[ localRow * columns() + localCol ] = value;
    }

    inline const RangeFieldType get ( const int localRow, 
                                      const int localCol ) const
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );
      return fields_[ localRow * columns() + localCol ];
    }

    /** \copydoc Dune::LocalMatrixInterface::clear */
    inline void clear ()
    {
      fields_.assign( 0 );
    }
  };
  
}

#endif
