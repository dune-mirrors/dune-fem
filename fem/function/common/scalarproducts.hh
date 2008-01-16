#ifndef DUNE_FEM_SCALARPRODURCTS_HH
#define DUNE_FEM_SCALARPRODURCTS_HH

//- system includes 
#include <iostream>
#include <map> 
#include <vector>

//- Dune includes  
#include <dune/common/mpihelper.hh>
#include <dune/grid/common/datahandleif.hh>

//- Dune-fem includes 
#include <dune/fem/space/common/commoperations.hh>
#include <dune/fem/space/common/singletonlist.hh>
#include <dune/fem/space/common/arrays.hh>
#include <dune/fem/space/common/gridpartutility.hh>

// in case of ISTL found include some more headers 
#if HAVE_DUNE_ISTL
#include <dune/istl/scalarproducts.hh>
#include <dune/fem/function/blockvectorfunction.hh>
#endif

namespace Dune
{
  
/** @addtogroup Communication Communication 
    @{
**/
  
  template< class SpaceImp, class MapperImp >
  class SlaveDofs 
  {
  private:
    // index map for send and receive data 
    class CommunicationIndexMap;

    template< class IndexMap, class Space >
    class LinkBuilder;

  public:
    //! type of discrete function space 
    typedef SpaceImp SpaceType; 
    //! type of grid part 
    typedef typename SpaceType :: GridPartType GridPartType;

    //! type of used mapper 
    typedef MapperImp MapperType;

  protected:
    const SpaceType &space_;
    const GridPartType &gridPart_;
    const MapperType &mapper_;

    const int myRank_;
    const int mySize_; 
    
    // type of communication indices 
    typedef CommunicationIndexMap IndexMapType;
    IndexMapType slaveDofs_;

    //! know grid sequence number 
    int sequence_; 
    
  public:
    //! constructor taking space 
    inline SlaveDofs ( const SpaceType &space,
                       const MapperType &mapper )
    : space_( space ),
      gridPart_( space_.gridPart() ),
      mapper_( mapper )
      myRank_( gridPart_.grid().comm().rank() ),
      mySize_( gridPart_.grid().comm().size() ),
      slaveDofs_(),
      sequence_( -1 )
    {}

  private:
    // prohibit copying
    SlaveDofs ( const SlaveDofs & );

  public:
    //! destrcutor removeing mpAccess 
    ~SlaveDofs()
    {}

    const int operator [] ( const int index ) const
    {
      return slaveDofs_[ index ];
    }

    const int size () const
    {
      return slaveDofs_.size();
    }

  public:
    inline void insert( const std :: vector< int > &indices )
    {
      slaveDofs_.insert( indices );
    }
    
    inline void initialize ()
    {
      sequence_ = -1;
      slaveDofs_.clear();
    }
    
    inline void finalize ()
    {
      // sort number for cache efficiency 
      slaveDofs_.sort();

      // store actual sequence number 
      sequence_ = space_.sequence();
    }
    
    //! check if grid has changed and rebuild cache if necessary 
    inline void rebuild () 
    {
      // check whether grid has changed. 
      if( sequence_ != space_.sequence() )
      {
        initialize();
        buildMaps();
        finalize();
      }
    }
    
  protected:  
    // build linkage and index maps 
    inline void buildMaps ()
    {
      if( !space_.continuous() )
        buildDiscontinuousMaps();
      else
        buildCommunicatedMaps();
    }
    
    // for discontinuous spaces we don't have to communicate
    inline void buildDiscontinuousMaps ();

    inline void buildCommunicatedMaps ();
  };


  
  template< class SpaceImp, class MapperImp >
  inline void SlaveDofs :: buildDiscontinuousMaps ()
  {
    typedef typename GridPartNewPartitionType< GridPartType, All_Partition >
      :: NewGridPartType NewGridPartType;
    typedef typename NewGridPartType :: template Codim<0> :: IteratorType IteratorType;

    NewGridPartType gridPart( const_cast< GridPartType & >( gridPart_ ).grid() );

    IteratorType endit = gridPart.template end<0>();
    for( IteratorType it = gridPart.template begin<0>(); it != endit; ++it )
    {
      typedef typename GridPartType :: GridType :: template Codim< 0 > :: Entity
        EntityType;
      
      const EntityType &entity = *it;
      if( entity.partitionType() != InteriorEntity ) 
      {
        // build local mapping 
        const int numDofs = mapper_.numEntityDofs( entity );
        std::vector< int > indices( numDofs );

        // copy numDofs 
        for( int i = 0; i < numDofs; ++i )
          indices[ i ] = mapper_.mapEntityDofToGlobal( entity, i );

        slaveDofs_.insert( indices ); 
      }
    }

    {
      // insert overall size at the end
      std :: vector< int > indices( 1, mapper_.size() );
      slaveDofs_.insert( indices );
    }

    //slaveDofs_.print(std::cout,myRank_);
  }



  template< class SpaceImp, class MapperImp >
  inline void SlaveDofs :: buildCommunicatedMaps ()
  {
    typedef LinkBuilder< IndexMapType, SpaceType > LinkBuilderHandleType; 
    LinkBuilderHandleType handle( slaveDofs_, space_ );

    gridPart_.communicate
      ( handle, InteriorBorder_All_Interface, ForwardCommunication );
    
    /*
    if( gridPart_.grid().overlapSize( 0 ) > 0 )
    {
      // case of YaspGrid, where overlap is treated as ghost
      gridPart_.communicate
        ( handle, InteriorBorder_All_Interface, ForwardCommunication );
      gridPart_.communicate
        ( handle, InteriorBorder_All_Interface, BackwardCommunication );
    }
    else 
    {
      // case of ALUGrid, where we have only the interior ghost situation 
      gridPart_.communicate( handle, All_All_Interface , ForwardCommunication );
    }
    */

    // insert overall size at the end
    std :: vector< int > indices( 1, mapper_.size() );
    slaveDofs_.insert( indices );
  }



  class SlaveDofs :: CommunicationIndexMap
  {
  protected:
    MutableArray< int > index_;

  public:
    //! constructor creating empty map
    CommunicationIndexMap()
    : index_( 0 )
    {
      index_.setMemoryFactor( 1.1 );
    }

  private:
    // prohibit copying
    CommunicationIndexMap( const CommunicationIndexMap& );

  public:
    //! reserve memory 
    void reserve( int size ) 
    {
      // resize array, memory factor will be used 
      index_.resize( size );
    }

    //! clear index map 
    void clear() 
    {
      // resize 0 will free memory 
      index_.resize( 0 );
    }

    //! append index vector with idx 
    void insert( const std :: vector< int > &idx )
    {
      const int size = idx.size();
      int count = index_.size();
      
      // reserve memory 
      reserve( count + size );
      assert( index_.size() == (count + size) );

      // copy indices to index vector 
      for( int i = 0; i < size; ++i, ++count )
      { 
        assert( idx[ i ] >= 0 );
        index_[ count ] = idx[ i ];
      }
    }

    //! return index map for entry i
    const int operator [] ( const int i ) const
    {
      assert( (i >= 0) && (i < size()) );
      return index_[ i ];
    }

    //! return size of map
    int size () const
    {
      return index_.size();
    }

    //! print  map for debugging only 
    void print( std :: ostream &s, int rank ) const
    {
      const int size = index_.size();
      s << "Start print: size = " << size << std :: endl;
      for( int i = 0; i < size; ++i )
        s << rank << " idx[ " << i << " ] = " << index_[ i ] << std :: endl;
      s << "End of Array" << std :: endl;
    }

    void sort() 
    {
      std :: sort( index_.begin(), index_.end() );
    }
  };



  template< class IndexMap, class Space >
  class SlaveDofs :: LinkBuilder
  : public CommDataHandleIF< LinkBuilder< IndexMap, Space >, int >
  {
  public:
    typedef IndexMap IndexMapType;
    typedef Space SpaceType;

    typedef typename SpaceType :: MapperType MapperType;

#if 0
  private:
    template< int dummy, int codim >
    struct CheckInterior
    {
      inline static bool check( const PartitionType p )
      {
        DUNE_THROW( NotImplemented, "Method not implemented!" );
        return true;
      }
    };

    // codim 0 specialization 
    template< int dummy >
    struct CheckInterior< dummy, 0 >
    {
      inline static bool check( const PartitionType p )
      {
        return (p == InteriorEntity);
      }
    };
#endif
    
  public:
    typedef int DataType;

    const int myRank_;
    const int mySize_;
    
    IndexMapType &slaveDofs_;

    const SpaceType &space_;
    const MapperType &mapper_;

  public:
    LinkBuilder( IndexMapType &slaveDofs,
                 const SpaceType &space )
    : myRank_( space.grid().comm().rank() ),
      mySize_( space.grid().comm().size() ),
      slaveDofs_( slaveDofs ),
      space_( space )
      mapper_( space.mapper() )
    {}

    bool contains ( int dim, int codim ) const
    {
      return space_.contains( codim );
    }

    bool fixedsize ( int dim, int codim ) const
    {
      return false;
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    inline void gather ( MessageBuffer &buffer,
                         const Entity &entity ) const
    {
      PartitionType ptype = entity.partitionType();

      if( (ptype == InteriorEntity) || (ptype == BorderEntity) )
        buff.write( myRank_ );
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    inline void scatter ( MessageBuffer &buffer,
                          const EntityType &entity,
                          size_t n )
    {
      PartitionType ptype = entity.partitionType();
      
      int minRank
        = ((ptype == InteriorEntity) || (ptype == BorderEntity) ? myRank_ : -1);
      for( size_t i = 0; i < n; ++i )
      {
        int rank;
        buff.read( rank );
        minRank = (rank < minRank ? rank : minRank);
      }

      // minimal rank means master
      assert( minRank != -1 );
      if( minRank != myRank_ )
      {
        // build local mapping 
        const int numDofs = mapper_.numEntityDofs( entity );
        std :: vector< int > indices( numDofs );

        // copy numDofs 
        for( int i = 0; i < numDofs; ++i )
          indices[ i ] = mapper_.mapEntityDofToGlobal( entity, i );
        
        // insert slave Dofs 
        slaveDofs_.insert( indices );
      }
    }

    //! return local dof size to be communicated 
    template< class Entity >
    size_t size ( const Entity &entity ) const
    {
      PartitionType ptype = entity.partitionType();

      return ((ptype == InteriorEntity) || (ptype == BorderEntity) : 1 : 0);
    }
  };

  

  //! Key for CommManager singleton list
  template< class SpaceImp, class MapperImp >
  class SlaveDofsSingletonKey
  {
  public:
    typedef SpaceImp SpaceType;
    typedef MapperImp MapperType;
    
  protected:
    const SpaceType &space_;
    const MapperType *const mapper_;

  public:
    //! constructor taking space 
    inline SlaveDofsSingletonKey ( const SpaceType &space, 
                                   const MapperType &mapper )
    : space_( space ),
      mapper_( &mapper )
    {}

    //! copy constructor  
    inline SlaveDofsSingletonKey ( const SlaveDofsSingletonKey &other )
    : space_( other.space_ ),
      mapper_( other.mapper_ )
    {}
    
    //! returns true if indexSet pointer and numDofs are equal 
    inline bool operator== ( const SlaveDofsSingletonKey &other ) const
    {
      return (space_ == other.space_) && (mapper_ == other.mapper_);
    }

    //! return reference to index set 
    const SpaceType &space () const
    {
      return space_;
    }

    //! return reference to index set 
    const MapperType &mapper () const
    {
      return *mapper_;
    }
  };



  //! Factory class for SingletonList to tell how objects are created and
  //! how compared.
  template< class KeyImp, class ObjectImp >
  class SlaveDofsFactory
  {
  public:
    //! create new communiaction manager   
    static ObjectImp *createObject( const KeyImp &key )
    {
      return new ObjectImp( key.space(), key.mapper() );
    }
    
    //! delete comm manager  
    static void deleteObject( ObjectImp *obj )
    {
      delete obj; 
    }
  };



  //! Proxy class to evaluate ScalarProduct 
  //! holding SlaveDofs which is singleton per space and mapper 
  template< class DiscreteFunction >
  class ParallelScalarProduct 
  {
  public:
    typedef DiscreteFunction DiscreteFunctionType;

  private:
    typedef ParallelScalarProduct< DiscreteFunctionType > ThisType;

  public:
    //! type of the discrete function space
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of range field 
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;

    //! type of used mapper 
    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    
    // type of communication manager object which does communication
    typedef SlaveDofs< DiscreteFunctionSpaceType, MapperType > SlaveDofsType;

    typedef SlaveDofsSingletonKey< DiscreteFunctionSpaceType, MapperType > KeyType;
    typedef SlaveDofsFactory< KeyType, SlaveDofsType > FactoryType;
    typedef SingletonList< KeyType, SlaveDofsType, FactoryType >
      SlaveDofsProviderType;

  protected:
    const DiscreteFunctionSpaceType &space_; 

    // is singleton per space 
    SlaveDofsType *const slaveDofs_;

  public:  
    //! constructor taking space 
    inline ParallelScalarProduct ( const DiscreteFunctionSpaceType &space )
    : space_( space ),
      slaveDofs_( getSlaveDofs( space_ ) )
    {
    }

  protected:
    inline static SlaveDofsType *getSlaveDofs ( DiscreteFunctionSpaceType &space )
    {
      KeyType key( space, space.mapper() );
      return &(SlaveDofsProviderType :: getObject( key ));
    }

  private:
    // prohibit copying
    ParallelScalarProduct( const ThisType & );

  public:
    //! remove object comm
    inline ~ParallelScalarProduct ()
    {
      SlaveDofsProviderType :: removeObject( *slaveDofs_ );
    }

    inline RangeFieldType scalarProductDofs ( const DiscreteFunctionType &x,
                                              const DiscreteFunctionType &y ) const
    {
      // rebuild slave dofs if grid was changed
      slaveDofs_->rebuild();
      const int numSlaves = slaveDofs_->size();

      RangeFieldType scp = 0;
      
      int i = 0;
      for( int slave = 0; slave < numSlaves; ++slave )
      {
        const int nextSlave = slaveDofs_[ slave ];
        for(; i < nextSlave; ++i )
          scp += x.dof( i ) * y.dof( i );
        // skip the slave dof
        ++i;
      }

      scp = space_.grid().comm().sum( scp );
      return scp;
    }
  };



#if HAVE_DUNE_ISTL
  //! Proxy class to evaluate ScalarProduct 
  //! holding SlaveDofs which is singleton per space and mapper 
  template <class DiscreteFunctionSpaceImp> 
  class ParallelScalarProduct<BlockVectorDiscreteFunction<DiscreteFunctionSpaceImp> > 
  : public ScalarProduct<typename BlockVectorDiscreteFunction<DiscreteFunctionSpaceImp> :: DofStorageType >
  {
    template<class SlaveDofsImp>
    class SlaveDofsProxy
    {
    protected:  
      SlaveDofsImp& slaveDofs_;
      SlaveDofsProxy(const SlaveDofsProxy& org); 
    public:
      SlaveDofsProxy(SlaveDofsImp& sd) : slaveDofs_(sd) 
      {
        slaveDofs_.initialize();
      }
      ~SlaveDofsProxy() 
      { 
        slaveDofs_.finalize(); 
      }

      void insert(const std::vector<int> & indices)
      {
        slaveDofs_.insert( indices );
      }
    };

    //! discrete function type 
    typedef BlockVectorDiscreteFunction<DiscreteFunctionSpaceImp> DiscreteFunctionType;
    //! type of this class 
    typedef ParallelScalarProduct<DiscreteFunctionType> ThisType;
    //! type of BlockVector 
    typedef typename DiscreteFunctionType :: DofStorageType BlockVectorType;
    //! type of discrete function space 
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of range field 
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;

    //! type of used mapper 
    typedef typename DiscreteFunctionSpaceType :: BlockMapperType MapperType;
    
  public:  
    // type of communication manager object which does communication 
    typedef SlaveDofs<DiscreteFunctionSpaceType,MapperType> SlaveDofsType;
  private:  

    typedef SlaveDofsSingletonKey<DiscreteFunctionSpaceType,MapperType> KeyType;
    typedef SlaveDofsFactory<KeyType, SlaveDofsType> FactoryType;

    typedef SingletonList< KeyType , SlaveDofsType , FactoryType > SlaveDofsProviderType;
  public:
    typedef SlaveDofsProxy<SlaveDofsType> BuildProxyType;
    //! export types
    typedef BlockVectorType domain_type;
    typedef typename BlockVectorType :: block_type :: field_type field_type;

    //! define the category
    enum { category=SolverCategory::sequential };

  protected:
    const DiscreteFunctionSpaceType & space_; 
    const KeyType key_;

    // is singleton per space 
    mutable SlaveDofsType & slaveDofs_;

    ParallelScalarProduct(const ThisType& org);
  public:  
    //! constructor taking space 
    ParallelScalarProduct(const DiscreteFunctionSpaceType& space) 
      : space_(space)
      , key_(space_,space_.blockMapper())
      , slaveDofs_(SlaveDofsProviderType::getObject(key_)) 
    {
    }

    //! remove object comm
    ~ParallelScalarProduct() 
    {
      SlaveDofsProviderType::removeObject(slaveDofs_);
    }

    std::auto_ptr<BuildProxyType> buildProxy() { return std::auto_ptr<BuildProxyType> (new BuildProxyType(slaveDofs_)); }

    const SlaveDofsType& slaveDofs() const { return slaveDofs_; }

    /*! \brief Dot product of two discrete functions. 
      It is assumed that the vectors are consistent on the interior+border
      partition.
     */
    RangeFieldType scalarProductDofs(const DiscreteFunctionType& x,
                                     const DiscreteFunctionType& y) const 
    {
      return scalarProductDofs(x.blockVector(),y.blockVector());
    }

    /*! \brief Dot product of two vectors. 
      It is assumed that the vectors are consistent on the interior+border
      partition.
     */
    virtual field_type dot (const BlockVectorType& x, 
                            const BlockVectorType& y)
    {
      return const_cast<ThisType&> (*this).scalarProductDofs(x,y);
    }

    /*! \brief Norm of a right-hand side vector. 
      The vector must be consistent on the interior+border partition
     */
    virtual double norm (const BlockVectorType& x)
    {
      return std::sqrt( const_cast<ThisType&> (*this).scalarProductDofs(x,x) );
    }
  protected:    
    /*! \brief Dot product of two block vectors. 
      It is assumed that the vectors are consistent on the interior+border
      partition.
     */
    RangeFieldType scalarProductDofs(const BlockVectorType& x,
                                     const BlockVectorType& y) const 
    {
#if HAVE_MPI
      // rebuild slave dofs if grid was changed  
      slaveDofs_.rebuild();

      RangeFieldType scp = 0;
      int i = 0;
      const int slaveSize = slaveDofs_.size();
      for(int slave = 0; slave<slaveSize; ++slave)
      {
        const int nextSlave = slaveDofs_[slave];
        for(; i<nextSlave; ++i) 
        {
          scp += x[i] * y[i];
        }
        // set i to next valid value 
        i = nextSlave + 1;
      }
      scp = space_.grid().comm().sum( scp );
      return scp;
#else 
      // return build-in scalar product 
      return x * y;
#endif
    }
  };
#endif

  //@}
} // end namespace Dune 
#endif
