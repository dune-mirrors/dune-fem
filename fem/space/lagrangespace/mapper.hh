#ifndef DUNE_LAGRANGESPACE_MAPPER_HH
#define DUNE_LAGRANGESPACE_MAPPER_HH

//- Dune includes 
#include <dune/common/geometrytype.hh>
#include <dune/common/exceptions.hh>

//- Dune-Fem includes 
#include <dune/fem/misc/codimmap.hh>
#include <dune/fem/space/common/dofmanager.hh>

//- local includes 
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

    typedef LagrangePointSet< GridPartType, polynomialOrder >
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
    , maxDofs_ (0)
    {
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
    virtual ~LagrangeMapper () {}

    //! get size (i.e. size of DoF vector)
    int size () const
    {
      return DimRange * indexSet_.size( dimension );
    }

    //! map a local DoF number to a global one
    template< class EntityType >
    int mapToGlobal ( const EntityType &entity, const int local ) const
    {
      const int coordinate = local % DimRange;
      const int localDof = local / DimRange;
      const int globalDof
        = indexSet_.template subIndex< dimension >( entity, localDof );
      return DimRange * globalDof + coordinate;
    }

    //! vector index of the hole
    int oldIndex ( int hole, int ) const
    {
      const int coordinate = hole % DimRange;
      const int setHole = hole / DimRange;
      const int setIndex = indexSet_.oldIndex( setHole, dimension );
      return setIndex * DimRange + coordinate;
    }

    //! vector index of data to be copied to the hole
    int newIndex ( int hole , int ) const
    {
      const int coordinate = hole % DimRange;
      const int setHole = hole / DimRange;
      const int setIndex = indexSet_.newIndex( setHole, dimension );
      return setIndex * DimRange + coordinate;
    }

    //! number of holes in DoF vector
    int numberOfHoles (int) const
    {
      return DimRange * indexSet_.numberOfHoles( dimension );
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

    typedef LagrangePointSet< GridPartType, polynomialOrder >
      LagrangePointSetType;
    typedef std :: map< const GeometryType, const LagrangePointSetType* >
      LagrangePointSetMapType;

    typedef typename LagrangePointSetType :: DofInfo DofInfo;

    typedef DofManager<GridType> DofManagerType;
    typedef DofManagerFactory<DofManagerType> DMFactoryType;

  private:
    template< unsigned int codim >
    class IndexSetCodimCallImp
    : public IndexSetCodimCall< GridPartType, codim >
    {
    };
    
    typedef CodimMap< dimension+1, IndexSetCodimCallImp >
      IndexSetCodimCallMapType;

  private:
    // reference to dof manager needed for debug issues 
    const DofManagerType& dm_;

    const IndexSetType &indexSet_;
    const IndexSetCodimCallMapType indexSetCodimCall_;
    
    LagrangePointSetMapType &lagrangePointSet_;

    unsigned int maxDofs_[ dimension+1 ];
    mutable unsigned int offset_[ dimension+1 ];
    mutable unsigned int oldOffSet_[ dimension+1 ];
    mutable unsigned int size_;
    unsigned int numDofs_;

    // for debugging only 
    mutable int sequence_;
  public:
    LagrangeMapper ( const GridPartType &gridPart,
                     LagrangePointSetMapType &lagrangePointSet )
    : dm_( DMFactoryType :: getDofManager(gridPart.grid()) )
    , indexSet_( gridPart.indexSet() )
    , lagrangePointSet_( lagrangePointSet )
    , sequence_( dm_.sequence() )
    {
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
        oldOffSet_[ codim ] = size_;
        size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
      }

      numDofs_ = 0;
      for(int codim = 0; codim <= dimension; ++codim )
      {
        numDofs_ += maxDofs_[codim];
      }
    }
    
    //! destructor 
    virtual ~LagrangeMapper ()
    {
    }

    //! return overall number of degrees of freedom 
    int size () const
    {
      return DimRange * size_;
    }

    //! map local dof of entity to global dof number 
    template< class EntityType >
    int mapToGlobal ( const EntityType &entity, const int local ) const
    {
      const int coordinate = local % DimRange;
      const int localDof = local / DimRange;
      
      // unsigned int codim, subEntity;
      const LagrangePointSetType *set
        = lagrangePointSet_[ entity.geometry().type() ];
      const DofInfo& dofInfo = set->dofInfo( localDof );
      
      const int subIndex
        = indexSetCodimCall_[ dofInfo.codim ]
            .subIndex( indexSet_, entity, dofInfo.subEntity );
      return DimRange * (offset_[ dofInfo.codim ] + subIndex) + coordinate;
    }

    //! return old dof number for given number of hole and codim (=block) 
    int oldIndex ( const int num, const int codim ) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / DimRange);
      // local number of dof is local 
      const int local = (num % DimRange);
      // codim to be revised 
      return DimRange * (oldOffSet_[codim] + indexSet_.oldIndex(newn,codim)) + local;
    }

    //! return old new number for given number of hole and codima (=block)
    int newIndex ( const int num , const int codim) const
    {
      // corresponding number of set is newn 
      const int newn  = static_cast<int> (num / DimRange);
      // local number of dof is local 
      const int local = (num % DimRange);
      // codim to be revised 
      return DimRange * (offset_[codim] + indexSet_.newIndex(newn,codim)) + local;
    }

    //! return number of holes for given codim 
    int numberOfHoles ( const int codim ) const
    {
      return (maxDofs_[ codim ] > 0) ? 
        (DimRange * indexSet_.numberOfHoles( codim )) : 0;
    }

    //! update mapper, i.e. calculate new insertion points 
    void update()
    {
      // assure that update is only called once per 
      // dof manager resize 
      if( sequence_ != dm_.sequence() )
      {
        // calculate new size 
        size_ = 0;
        for( int codim = 0; codim <= dimension; ++codim )
        {
          oldOffSet_[ codim ] = offset_[codim];
          offset_[ codim ] = size_;
          size_ += indexSet_.size( codim ) * maxDofs_[ codim ];
        }
        sequence_ = dm_.sequence();
      }
      else 
      {
        DUNE_THROW(InvalidStateException,"update of mapper should only be called once per sequence");
      }
    }

    //! return number of supported codims 
    int numBlocks () const 
    {
      return dimension + 1;
    }

    //! return old offsets for block 
    int oldOffSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return DimRange * oldOffSet_[ block ];
    }

    //! return current offsets 
    int offSet ( const int block ) const
    {
      assert( (block >= 0) && (block < numBlocks()) );
      return DimRange * offset_[ block ];
    }

    //! return number of dofs per element 
    int numDofs () const
    {
      return numDofs_;
    }

    //! calculate new size without changing object 
    int newSize () const
    {
      int newSize = 0;
      for( int codim = 0; codim <= dimension; ++codim )
      {
        newSize += indexSet_.size( codim ) * maxDofs_[ codim ];
      }
      return DimRange * newSize;
    }

    //! returns true if compression is needed 
    bool needsCompress () const
    {
      return indexSet_.needsCompress();
    }
  };
 
} // end namespace Dune 
#endif
