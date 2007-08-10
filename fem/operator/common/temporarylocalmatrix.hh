#ifndef DUNE_FEM_TEMPORARYLOCALMATRIX_HH
#define DUNE_FEM_TEMPORARYLOCALMATRIX_HH

#include <dune/fem/misc/array.hh>
#include <dune/fem/operator/common/localmatrix.hh>

namespace Dune
{

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

    template< class DomainEntityType, class RangeEntityType >
    inline void init ( const DomainEntityType &domainEntity,
                       const RangeEntityType &rangeEntity )
    {
      BaseType :: init( domainEntity, rangeEntity );
      fields_.resize( rows() * columns() );
    }

    inline void add( const int localRow, 
                     const int localCol, 
                     const RangeFieldType &value )
    {
      assert( (localRow >= 0) && (localRow < rows()) );
      assert( (localCol >= 0) && (localCol < columns()) );
      fields_[ localRow * columns() + localCol ] += value;
    }

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

    inline void clear ()
    {
      fields_.assign( 0 );
    }
  };
  
}

#endif
