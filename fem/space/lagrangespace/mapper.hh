#ifndef DUNE_LAGRANGESPACE_MAPPER_HH
#define DUNE_LAGRANGESPACE_MAPPER_HH

#include <dune/common/geometrytype.hh>

#include "lagrangepoints.hh"
#include "indexsetcodimcall.hh"

namespace Dune
{
  
  template< class GridPartImp, unsigned int polOrder, unsigned int dimrange >
  class LagrangeMapper;

  template< class GridPartImp, unsigned int dimrange >
  class LagrangeMapper< GridPartImp, 1, dimrange >
  : public DofMapperDefault< LagrangeMapper< GridPartImp, 1, dimrange > >
  {
  public:
    typedef GridPartImp GridPartType;
    
    typedef typename GridPartType :: GridType GridType;

    typedef typename GridType :: ctype FieldType;

    enum { dimension = GridType :: dimension };

    enum { polynomialOrder = 1 };

    enum { DimRange = dimrange };

  private:
    typedef LagrangeMapper< GridPartType, polynomialOrder, DimRange > ThisType;
    typedef DofMapperDefault< ThisType > BaseType;

  public:
    typedef typename GridPartType :: IndexSetType IndexSetType;

    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      LagrangePointSetType;
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

  private:
    const IndexSetType &indexSet_;
    unsigned int maxDofs_;

  public:
    //! constructor
    LagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapType &lagrangePointSet )
    : indexSet_( gridPart.indexSet() )
    {
      maxDofs_ = 0;
        
      typedef typename LagrangePointSetMapType :: iterator IteratorType;
      IteratorType end = lagrangePointSet.end();
      for( IteratorType it = lagrangePointSet.begin(); it != end; ++it ) {
        const LagrangePointSetType *set = (*it).second;
        if( set == NULL )
          continue;
        
        const unsigned int setDofs = set->maxDofs( dimension );
        maxDofs_ = (maxDofs_ >= setDofs) ? maxDofs_ : setDofs;
      }
    }
   
    //! destructor
    virtual ~LagrangeMapper ()
    {
    }

    //! get size (i.e. size of DoF vector)
    int size () const
    {
      return DimRange * indexSet_.size( dimension );
    }

    //! map a local DoF number to a global one
    template< class EntityType >
    int mapToGlobal ( EntityType &entity, int local ) const
    {
      const int coordinate = local % DimRange;
      const int localDof = local / DimRange;
      const int globalDof
        = indexSet_.template index< dimension >( entity, localDof );
      return DimRange * globalDof + coordinate;
    }

    //! vector index of the hole
    int oldIndex ( int hole ) const
    {
      const int coordinate = hole % DimRange;
      const int setHole = hole / DimRange;
      const int setIndex = indexSet_.oldIndex( setHole, dimension );
      return setIndex * DimRange + coordinate;
    }

    //! vector index of data to be copied to the hole
    int newIndex ( int hole ) const
    {
      const int coordinate = hole % DimRange;
      const int setHole = hole / DimRange;
      const int setIndex = indexSet_.newIndex( setHole, dimension );
      return setIndex * DimRange + coordinate;
    }

    //! number of holes in DoF vector
    int numberOfHoles () const
    {
      return DimRange * indexSet_.numberOfHoles( dimension );
    }

    //! additional size needed during restriction
    int additionalSizeEstimate () const
    {
      return DimRange * indexSet_.additionalSizeEstimate();
    }

    //! get maximum number of DoFs per entity
    int numDofs () const
    {
      return DimRange * maxDofs_;
    }

    //! get size of DoF vector after adaption
    int newSize () const
    {
      return this->size();
    }

    //! should the DoF vector be compressed?
    bool needsCompress () const
    {
      return indexSet_.needsCompress();
    }
  };


  template< class GridPartImp, unsigned int dimrange >
  class LagrangeMapper< GridPartImp, 2, dimrange >
  : public DofMapperDefault< LagrangeMapper< GridPartImp, 2, dimrange > >
  {
  public:
    typedef GridPartImp GridPartType;
    
    typedef typename GridPartType :: GridType GridType;

    typedef typename GridType :: ctype FieldType;

    enum { dimension = GridType :: dimension };

    enum { polynomialOrder = 2 };

    enum { DimRange = dimrange };

  private:
    typedef LagrangeMapper< GridPartType, polynomialOrder, DimRange > ThisType;
    typedef DofMapperDefault< ThisType > BaseType;

  public:
    typedef typename GridPartType :: IndexSetType IndexSetType;

    typedef LagrangePointSetInterface< FieldType, dimension, polynomialOrder >
      LagrangePointSetType;
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    typedef IndexSetCodimCallFactory< GridPartType >
      IndexSetCodimCallFactoryType;
    typedef typename IndexSetCodimCallFactoryType :: IndexSetCodimCallType
      IndexSetCodimCallType;

  private:
    const IndexSetType &indexSet_;
    IndexSetCodimCallType *indexSetCodimCall_[ dimension+1 ];
    
    LagrangePointSetMapType &lagrangePointSet_;

    unsigned int maxDofs_[ dimension+1 ];
    unsigned int offset_[ dimension+1 ];
    unsigned int size_;

  public:
    LagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapType &lagrangePointSet )
    : indexSet_( gridPart.indexSet() ),
      lagrangePointSet_( lagrangePointSet )
    {
      IndexSetCodimCallFactoryType :: getAllCalls( indexSetCodimCall_ );
        
      for( int codim = 0; codim <= dimension; ++codim )
        maxDofs_[ codim ] = 0;
      
      typedef typename LagrangePointSetMapType :: iterator IteratorType;
      IteratorType end = lagrangePointSet_.end();
      for( IteratorType it = lagrangePointSet_.begin(); it != end; ++it ) {
        const LagrangePointSetType *set = (*it).second;
        if( set == NULL )
          continue;
        
        for( int codim = 0; codim <= dimension; ++codim ) {
          const unsigned int setDofs = set->maxDofs( codim );
          unsigned int &maxDofs = maxDofs_[ codim ];
          maxDofs = (maxDofs >= setDofs) ? maxDofs : setDofs;
        }
      }

      size_ = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        offset_[ codim ] = size_;
        size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
      }
    }
    
    virtual ~LagrangeMapper ()
    {
    }

    int size () const
    {
      return dimrange * size_;
    }

    template< class EntityType >
    int mapToGlobal ( EntityType &entity, int local ) const
    {
      const int coordinate = local % DimRange;
      const int localDof = local / DimRange;
        
      unsigned int codim, subEntity;
      const LagrangePointSetType *set
        = lagrangePointSet_[ entity.geometry().type() ];
      set->dofSubEntity( localDof, codim, subEntity );

      const int subIndex
        = indexSetCodimCall_[ codim ]->subIndex( indexSet_, entity, subEntity );
      return DimRange * (offset_[ codim ] + subIndex) + coordinate;
    }

    int oldIndex ( int elNum ) const
    {
      DUNE_THROW( NotImplemented, "LagrangeMapper not implemented yet." );
    }

    int newIndex ( int elNum ) const
    {
      DUNE_THROW( NotImplemented, "LagrangeMapper not implemented yet." );
    }

    int numberOfHoles () const
    {
      DUNE_THROW( NotImplemented, "LagrangeMapper not implemented yet." );
    }

    int additionalSizeEstimate () const
    {
      DUNE_THROW( NotImplemented, "LagrangeMapper not implemented yet." );
    }

    int numDofs () const
    {
      DUNE_THROW( NotImplemented, "LagrangeMapper not implemented yet." );
    }

    int newSize () const
    {
      DUNE_THROW( NotImplemented, "LagrangeMapper not implemented yet." );
    }

    bool needsCompress () const
    {
      DUNE_THROW( NotImplemented, "LagrangeMapper not implemented yet." );
    }
  };
 
}

#endif
