#ifndef DUNEFEM_SCALARPRODURCTS_HH
#define DUNEFEM_SCALARPRODURCTS_HH

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

namespace Dune { 
  
/** @addtogroup Communication Communication 
    @{
**/
  
  template <class SpaceImp, class MapperImp> 
  class SlaveDofs 
  {
    // index map for send and receive data 
    class CommunicationIndexMap
    {
      MutableArray<int> index_;
      CommunicationIndexMap(const CommunicationIndexMap&);
    public:
      //! constructor creating empty map
      CommunicationIndexMap() : index_(0) 
      {
        index_.setMemoryFactor(1.1);
      }

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
      void insert( const std::vector<int> & idx )
      {
        const int size = idx.size();
        int count = index_.size();
        // reserve memory 
        reserve( count + size );
        assert( index_.size() == (count+size));
        // copy indices to index vector 
        for(int i=0; i<size; ++i, ++count) 
        { 
          assert( idx[i] >= 0 );
          index_[count] = idx[i]; 
        }
      }

      //! return index map for entry i
      const int operator [] (const int i) const 
      {
        assert( i >= 0 );
        assert( i < (int) index_.size());
        return index_[i];
      }

      //! return size of map
      int size () const { return index_.size(); }

      //! print  map for debugging only 
      void print(std::ostream & s, int rank) const 
      {
        const int size = index_.size();
        s << "Start print: size = " << size << std::endl;
        for(int i=0; i<size; ++i) 
        {
          s<< rank << " idx["<<i<<"] = " << index_[i] << std::endl;
        }
        s << "End of Array\n";
      }

      void sort() 
      {
        std::sort( index_.begin(), index_.end() );
      }
    };

    template <class IndexMapType, 
              class SpaceType> 
    class LinkBuilder
     : public CommDataHandleIF<
       LinkBuilder< IndexMapType , SpaceType > ,
       int >
    {
      template <int dummy, int codim> 
      struct CheckInterior 
      {
        inline static bool check(const PartitionType p) 
        {
          DUNE_THROW(NotImplemented,"Method not implemented!");
          return true; 
        }
      };

      // codim 0 specialization 
      template <int dummy> 
      struct CheckInterior<dummy,0>
      {
        inline static bool check(const PartitionType p) 
        {
          return (p == InteriorEntity); 
        }
      };
      
    public:
      typedef int DataType;

      const int myRank_;
      const int mySize_; 
      
      IndexMapType & slaveDofs_;

      const SpaceType & space_;
      typedef typename SpaceType :: MapperType MapperType; 
      const MapperType& mapper_; 

    public:
      LinkBuilder(IndexMapType& slaveDofs,
                  const SpaceType & space)
        : myRank_(space.grid().comm().rank()) 
        , mySize_(space.grid().comm().size())
        , slaveDofs_(slaveDofs)
        , space_(space)
        , mapper_(space.mapper())
      {
      }

      bool contains (int dim, int codim) const
      {
        return space_.contains(codim);
      }

      bool fixedsize (int dim, int codim) const
      {
        return false;
      }

      //! read buffer and apply operation 
      template<class MessageBufferImp, class EntityType>
      void gather (MessageBufferImp& buff, const EntityType& en) const
      {
        // send rank 
        buff.write( myRank_ );
      }

      //! read buffer and apply operation 
      template<class MessageBufferImp, class EntityType>
      void scatter (MessageBufferImp& buff, const EntityType& en, size_t n)
      {
        int rank; 
        int minRank = myRank_;
        for(size_t i=0; i<n; ++i) 
        {
          buff.read( rank );
          if( rank < minRank ) minRank = rank;
        }

        // minimal rank means master 
        if( minRank != myRank_ || en.partitionType() != InteriorEntity ) 
        {
          // build local mapping 
          const int numDofs = mapper_.numEntityDofs( en );
          std::vector<int> indices(numDofs);

          // copy numDofs 
          for(int i=0; i<numDofs; ++i) 
          {
            indices[i] = mapper_.mapEntityDofToGlobal( en , i ); 
          }
          
          // insert slave Dofs 
          slaveDofs_.insert( indices );
        }
      }

      //! return local dof size to be communicated 
      template<class EntityType>
      size_t size (const EntityType& en) const
      {
        // size of data 
        return 1; 
      }
    };

    //! type of discrete function space 
    typedef SpaceImp SpaceType; 
    //! type of grid part 
    typedef typename SpaceType :: GridPartType GridPartType; 

    //! type of used mapper 
    typedef MapperImp MapperType;

    const SpaceType & space_; 
    const GridPartType & gridPart_; 
    const MapperType& mapper_;

    const int myRank_;
    const int mySize_; 
    
    // type of communication indices 
    typedef CommunicationIndexMap IndexMapType;
    IndexMapType slaveDofs_;
    // type of IndexMapVector 

    //! know grid sequence number 
    int sequence_; 
    
    // do not copy this class 
    SlaveDofs(const SlaveDofs &);
  public:
    //! constructor taking space 
    SlaveDofs(const SpaceType & space, const MapperType& mapper)
      : space_(space) , gridPart_(space_.gridPart()) 
      , mapper_(mapper)
      , myRank_(gridPart_.grid().comm().rank())
      , mySize_(gridPart_.grid().comm().size())
      , slaveDofs_()
      , sequence_(-1)
    {
    }

    //! destrcutor removeing mpAccess 
    ~SlaveDofs()
    {
    }

    const int operator [] (const int index) const 
    {
      return slaveDofs_[index];
    }

    const int size () const { return slaveDofs_.size(); }

  public:
    // build linkage and index maps 
    void buildMaps() 
    {
      slaveDofs_.clear();

      if( ! space_.continuous() )
      {
        typedef typename GridPartNewPartitionType<GridPartType,All_Partition>:: NewGridPartType NewGridPartType;
        typedef typename NewGridPartType :: template Codim<0> :: IteratorType IteratorType;

        NewGridPartType gridPart( const_cast<GridPartType&> (gridPart_).grid() );

        IteratorType endit = gridPart. template end<0>();
        for(IteratorType it = gridPart. template begin<0>();
            it != endit; ++it)
        {
          typedef typename GridPartType :: GridType :: template Codim<0>::Entity EntityType;
          const EntityType& en = *it; 
          if(en.partitionType() != InteriorEntity ) 
          {
            // build local mapping 
            const int numDofs = mapper_.numEntityDofs( en );
            std::vector<int> indices(numDofs);

            // copy numDofs 
            for(int i=0; i<numDofs; ++i)
            {
              indices[i] = mapper_.mapEntityDofToGlobal( en , i );
            }

            slaveDofs_.insert( indices ); 
          }
        }

        // insert overall size 
        {
          std::vector<int> indices(1);
          indices[0] = mapper_.size();
          slaveDofs_.insert( indices );
        }

        slaveDofs_.print(std::cout,0);
      }
      else 
      {
        DUNE_THROW(NotImplemented,"Do some work here, Martin!");

        /*
        // type of data handler 
        typedef LinkBuilder< IndexMapType , SpaceType > LinkBuilderHandleType; 
        LinkBuilderHandleType handle( slaveDofs_ , space_ );

        // do communication to build up linkage 
        if( gridPart_.grid().overlapSize(0) > 0 ) 
        {
          // case of YaspGrid, where Overlap is treated as ghost 
          gridPart_.communicate( handle, InteriorBorder_All_Interface , ForwardCommunication );
          gridPart_.communicate( handle, InteriorBorder_All_Interface , BackwardCommunication );
        }
        else 
        {
          // case of ALUGrid, where we have only the interior ghost situation 
          gridPart_.communicate( handle, All_All_Interface , ForwardCommunication );
        }
        */
      }
    }

    //! check if grid has changed and rebuild cache if necessary 
    void rebuild() 
    {
      // check whether grid has changed. 
      if(sequence_ != space_.sequence()) 
      {
        buildMaps();

        // store actual sequence number 
        sequence_ = space_.sequence();
      }
    }
  };

  //! Key for CommManager singleton list 
  template <class SpaceImp, class MapperImp>
  class SlaveDofsSingletonKey
  {
    const SpaceImp& space_;
    const MapperImp& mapper_;
  public:
    //! constructor taking space 
    SlaveDofsSingletonKey(const SpaceImp & space, 
                          const MapperImp& mapper) 
      : space_(space) , mapper_(mapper) 
    {}

    //! copy constructor  
    SlaveDofsSingletonKey(const SlaveDofsSingletonKey & org) 
      : space_(org.space_), mapper_(org.mapper_)
    {}
    //! returns true if indexSet pointer and numDofs are equal 
    bool operator == (const SlaveDofsSingletonKey & otherKey) const
    {
      // mapper of space is singleton 
      return ( (&(space_.mapper()) == & (otherKey.space_.mapper()) ) && 
               (&mapper_ == &otherKey.mapper_) 
             );
    }

    //! return reference to index set 
    const SpaceImp & space() const { return space_; }

    //! return reference to index set 
    const MapperImp & mapper() const { return mapper_; }
  };

  //! Factory class for SingletonList to tell how objects are created and
  //! how compared.
  template <class KeyImp, class ObjectImp>
  class SlaveDofsFactory
  {
    public:
    //! create new communiaction manager   
    static ObjectImp * createObject( const KeyImp & key )
    {
      return new ObjectImp(key.space(),key.mapper());
    }
    
    //! delete comm manager  
    static void deleteObject( ObjectImp * obj )
    {
      delete obj; 
    }
  };

  //! Proxy class to evaluate ScalarProduct 
  //! holding SlaveDofs which is singleton per space and mapper 
  template <class DiscreteFunctionImp> 
  class ParallelScalarProduct 
  {
    typedef ParallelScalarProduct<DiscreteFunctionImp> ThisType;
    typedef DiscreteFunctionImp DiscreteFunctionType;
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    //! type of range field 
    typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;

    //! type of used mapper 
    typedef typename DiscreteFunctionSpaceType :: MapperType MapperType;
    
    // type of communication manager object which does communication 
    typedef SlaveDofs<DiscreteFunctionSpaceType,MapperType> SlaveDofsType;

    typedef SlaveDofsSingletonKey<DiscreteFunctionSpaceType,MapperType> KeyType;
    typedef SlaveDofsFactory<KeyType, SlaveDofsType> FactoryType;

    typedef SingletonList< KeyType , SlaveDofsType , FactoryType > SlaveDofsProviderType;

    const DiscreteFunctionSpaceType & space_; 

    const KeyType key_;

    // is singleton per space 
    mutable SlaveDofsType & slaveDofs_;

    ParallelScalarProduct(const ThisType& org);
  public:  
    //! constructor taking space 
    ParallelScalarProduct(const DiscreteFunctionSpaceType& space) 
      : space_(space)
      , key_(space_,space_.mapper())
      , slaveDofs_(SlaveDofsProviderType::getObject(key_)) 
    {
    }

    //! remove object comm
    ~ParallelScalarProduct() 
    {
      SlaveDofsProviderType::removeObject(slaveDofs_);
    }

    RangeFieldType scalarProductDofs(const DiscreteFunctionType& x,
                                     const DiscreteFunctionType& y) const 
    {
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
          scp += x.dof(i) * y.dof(i);
        }
        // set i to next valid value 
        i = nextSlave + 1;
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
    
    // type of communication manager object which does communication 
    typedef SlaveDofs<DiscreteFunctionSpaceType,MapperType> SlaveDofsType;

    typedef SlaveDofsSingletonKey<DiscreteFunctionSpaceType,MapperType> KeyType;
    typedef SlaveDofsFactory<KeyType, SlaveDofsType> FactoryType;

    typedef SingletonList< KeyType , SlaveDofsType , FactoryType > SlaveDofsProviderType;
  public:
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
//#if HAVE_MPI
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
//#else 
//      // return build-in scalar product 
//      return x * y;
//#endif
    }
  };
#endif

  //@}
} // end namespace Dune 
#endif
