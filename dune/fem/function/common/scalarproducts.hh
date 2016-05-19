#ifndef DUNE_FEM_SCALARPRODURCTS_HH
#define DUNE_FEM_SCALARPRODURCTS_HH

#include <iostream>
#include <memory>
#include <set>
#include <map>
#include <limits>
#include <algorithm>

#include <dune/common/exceptions.hh>
#include <dune/common/genericiterator.hh>
#include <dune/common/ftraits.hh>
#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#if HAVE_DUNE_ISTL
#include <dune/istl/scalarproducts.hh>
#endif

#include <dune/fem/storage/singletonlist.hh>
#include <dune/fem/misc/mpimanager.hh>
#include <dune/fem/space/common/commindexmap.hh>
#include <dune/fem/function/blockvectorfunction/declaration.hh>

namespace Dune
{

  namespace Fem
  {

  /** @addtogroup Communication Communication
      @{
  **/

    template< class Space, class Mapper >
    class SlaveDofs
    {
      typedef SlaveDofs< Space, Mapper > ThisType;

    public:
      class SingletonKey;

    private:
      class LinkBuilder;

    public:
      //! type of discrete function space
      typedef Space SpaceType;

      //! for convenience
      typedef SpaceType DiscreteFunctionSpaceType ;

      //! type of grid part
      typedef typename SpaceType :: GridPartType GridPartType;

      //! type of used mapper
      typedef Mapper MapperType;

    protected:
      typedef Fem :: CommunicationIndexMap IndexMapType;

      const SpaceType &space_;
      const GridPartType &gridPart_;
      const MapperType &mapper_;

      const int myRank_;
      const int mySize_;

      // type of communication indices
      IndexMapType slaves_;
      std::set<int> slaveSet_;

      //! know grid sequence number
      int sequence_;

    public:
      //! constructor taking space
      SlaveDofs ( const SingletonKey &key )
      : space_( key.space() ),
        gridPart_( space_.gridPart() ),
        mapper_( key.mapper() ),
        myRank_( gridPart_.comm().rank() ),
        mySize_( gridPart_.comm().size() ),
        slaves_(),
        sequence_( -1 )
      {}

      SlaveDofs ( const SlaveDofs& ) = delete;

      //! return dof number of salve with index
      int operator [] ( const int index ) const
      {
        return slaves_[ index ];
      }

      //! return number of slave dofs
      int size () const
      {
        return slaves_.size();
      }

      //! return true if index is contained, meaning is a slave dof
      bool isSlave( const int index ) const
      {
        typedef GenericIterator<const IndexMapType, const int> IteratorType;

        return std::binary_search(IteratorType(slaves_, 0),
                                  IteratorType(slaves_, size()-1),
                                  index);
      }

      //! insert index
      void insert( const int index )
      {
        slaveSet_.insert( index );
      }

      //! initialize
      void initialize ()
      {
        sequence_ = -1;
        slaveSet_.clear();
        slaves_.clear();
      }

      //! finalize
      void finalize ()
      {
        // insert slaves
        slaves_.set( slaveSet_ );

        // remove memory
        slaveSet_.clear();

        // store actual sequence number
        sequence_ = space_.sequence();
      }

      //! check if grid has changed and rebuild cache if necessary
      void rebuild ()
      {
        // check whether grid has changed.
        if( sequence_ != space_.sequence() )
        {
          initialize();
          buildMaps();
          finalize();
        }
      }

      //! return reference to discrete function space
      const SpaceType& space () const
      {
        return space_;
      }

    protected:
      void buildMaps ()
      {
        // build linkage and index maps
        if( !space_.continuous() )
          buildDiscontinuousMaps();
        else
          buildCommunicatedMaps();
      }

      void buildDiscontinuousMaps ()
      {
        // for discontinuous spaces we don't have to communicate
        const auto idxpitype = GridPartType :: indexSetPartitionType;
        const auto endit = gridPart_.template end< 0, idxpitype >();
        for( auto it = gridPart_.template begin< 0, idxpitype >(); it != endit; ++it )
        {
          const auto& entity = *it;
          if( entity.partitionType() != Dune::InteriorEntity )
            mapper_.mapEachEntityDof( entity, [this]( const int, const auto& value ){this->insert( value );} );
        }
        // insert overall size at the end
        insert( mapper_.size() );
      }

      void buildCommunicatedMaps ()
      {
        // we have to skip communication when parallel program is executed only on one processor
        // otherwise YaspGrid and Lagrange polorder=2 fails :(
        if( gridPart_.comm().size() > 1 )
        {
          try
          {
            LinkBuilder handle( *this, space_ , mapper_ );
            gridPart_.communicate( handle, GridPartType::indexSetInterfaceType, ForwardCommunication );
          }
          catch( const Exception &e )
          {
            std::cerr << e << std::endl << "Exception thrown in: " << __FILE__ << " line:" << __LINE__ << std::endl;
            abort();
          }
        }
        // insert overall size at the end
        insert( mapper_.size() );
      }
    };



    //! Key for CommManager singleton list
    template< class Space, class Mapper >
    class SlaveDofs< Space, Mapper > :: SingletonKey
    {
    public:
      typedef Space SpaceType;
      typedef Mapper MapperType;

    protected:
      const SpaceType &space_;
      const MapperType *const mapper_;

    public:
      //! constructor taking space
      SingletonKey ( const SpaceType &space, const MapperType &mapper )
      : space_( space ), mapper_( &mapper )
      {}

      //! copy constructor
      SingletonKey ( const SingletonKey &other )
      : space_( other.space_ ), mapper_( other.mapper_ )
      {}

      //! returns true if indexSet pointer and numDofs are equal
      bool operator== ( const SingletonKey &other ) const
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



    template< class Space, class Mapper >
    class SlaveDofs< Space,Mapper > :: LinkBuilder
    : public CommDataHandleIF< LinkBuilder, int >
    {
    public:
      typedef Space SpaceType;
      typedef Mapper MapperType;

      enum { nCodim = SpaceType :: GridType :: dimension + 1 };

      typedef int DataType;

      const int myRank_;
      const int mySize_;

      typedef SlaveDofs< Space,Mapper > IndexMapType;
      IndexMapType &slaves_;

      const SpaceType &space_;
      const MapperType &mapper_;

      LinkBuilder( IndexMapType &slaves, const SpaceType &space, const MapperType& mapper )
      : myRank_( space.gridPart().comm().rank() ), mySize_( space.gridPart().comm().size() ),
        slaves_( slaves ), space_( space ), mapper_( mapper )
      {}

      bool contains ( int dim, int codim ) const
      {
        return mapper_.contains( codim );
      }

      bool fixedsize ( int dim, int codim ) const
      {
        return false;
      }

      //! read buffer and apply operation
      template< class MessageBuffer, class Entity >
      void gather ( MessageBuffer &buffer, const Entity &entity ) const
      {
        // for sending ranks write rank
        if( sendRank( entity ) )
          buffer.write( myRank_ );
      }

      //! read buffer and apply operation
      //! scatter is called for one every entity
      //! several times depending on how much data
      //! was gathered
      template< class MessageBuffer, class EntityType >
      void scatter ( MessageBuffer &buffer, const EntityType &entity, size_t n )
      {
        if( n == 1 )
        {
          int rank;
          buffer.read( rank );

          assert( (rank >= 0) && (rank < mySize_) );

          // if entity in not interiorBorder insert anyway
          if ( rank < myRank_ || ! sendRank( entity ) )
            mapper_.mapEachEntityDof( entity, [this]( const int , const auto& value ){slaves_.insert( value );} );
        }
      }

      //! return local dof size to be communicated
      template< class Entity >
      size_t size ( const Entity &entity ) const
      {
        return (sendRank( entity )) ? 1 : 0;
      }

    protected:
      template <class Entity>
      bool sendRank(const Entity& entity) const
      {
        const PartitionType ptype = entity.partitionType();
        return (ptype == InteriorEntity) || (ptype == BorderEntity);
      }
    };

    //! Proxy class to evaluate ScalarProduct
    //! holding SlaveDofs which is singleton per space and mapper
    template< class DiscreteFunctionSpace >
    class SlaveDofsProvider
    {
    public:
      //! type of the discrete function space
      typedef DiscreteFunctionSpace DiscreteFunctionSpaceType;

    private:
      typedef SlaveDofsProvider< DiscreteFunctionSpaceType > ThisType;

    public:
      //! type of used mapper
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType MapperType;

      enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

      // type of communication manager object which does communication
      typedef SlaveDofs< DiscreteFunctionSpaceType, MapperType > SlaveDofsType;
      typedef typename SlaveDofsType :: SingletonKey SlaveDofsKeyType;

      typedef SingletonList< SlaveDofsKeyType, SlaveDofsType >
        SlaveDofsProviderType;

    protected:
      const DiscreteFunctionSpaceType &space_;

      // is singleton per space
      SlaveDofsType *const slaveDofs_;

    public:
      //! constructor taking space
      SlaveDofsProvider ( const DiscreteFunctionSpaceType &space )
      : space_( space ), slaveDofs_( getSlaveDofs( space_ ) )
      {}

      SlaveDofsProvider( const ThisType& ) = delete;

      //! return discrete function space
      const DiscreteFunctionSpaceType& space() const
      {
        return space_;
      }

      //! remove object comm
      ~SlaveDofsProvider ()
      {
        SlaveDofsProviderType :: removeObject( *slaveDofs_ );
      }

      SlaveDofsType &slaveDofs () const
      {
        // rebuild slave dofs if grid was changed
        slaveDofs_->rebuild();
        return *slaveDofs_;
      }

    protected:
      static SlaveDofsType *getSlaveDofs ( const DiscreteFunctionSpaceType &space )
      {
        SlaveDofsKeyType key( space, space.blockMapper() );
        return &(SlaveDofsProviderType :: getObject( key ));
      }
    };

#if HAVE_MPI
    //! Proxy class to evaluate ScalarProduct
    //! holding SlaveDofs which is singleton per space and mapper
    template< class DiscreteFunction >
    class ParallelScalarProduct
      : public SlaveDofsProvider< typename DiscreteFunction :: DiscreteFunctionSpaceType >
    {
    public:
      typedef DiscreteFunction DiscreteFunctionType;

      //! type of the discrete function space
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
        DiscreteFunctionSpaceType;

    private:
      typedef ParallelScalarProduct< DiscreteFunctionType > ThisType;
      typedef SlaveDofsProvider< DiscreteFunctionSpaceType > BaseType;

    public:
      //! type of range field
      typedef typename DiscreteFunctionSpaceType :: RangeFieldType  RangeFieldType;

      //! type of used mapper
      typedef typename DiscreteFunctionSpaceType :: BlockMapperType MapperType;

      enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };

      // type of communication manager object which does communication
      typedef SlaveDofs< DiscreteFunctionSpaceType, MapperType > SlaveDofsType;

      typedef typename DiscreteFunctionType :: ConstDofBlockPtrType
        ConstDofBlockPtrType;

      //! constructor taking space
      ParallelScalarProduct ( const DiscreteFunctionSpaceType &space )
      : BaseType( space )
      {}

      using BaseType :: space;

      //! evaluate scalar product and omit slave nodes
      template < class OtherDiscreteFunctionType >
      RangeFieldType scalarProductDofs ( const DiscreteFunctionType &x, const OtherDiscreteFunctionType &y ) const
      {
        auto &slaveDofs = this->slaveDofs();

        RangeFieldType scp = 0;

        const int numSlaves = slaveDofs.size();
        for( int slave = 0, i = 0 ; slave < numSlaves; ++slave )
        {
          const int nextSlave = slaveDofs[ slave ];
          for(; i < nextSlave; ++i )
            for( unsigned int j = 0; j < blockSize; ++j )
              scp += x.dofVector()[ i ][ j ] * y.dofVector()[ i ][ j ];

          // skip the slave dof
          ++i;
        }

        // do global sum
        scp = space().gridPart().comm().sum( scp );
        return scp;
      }
    };
#else
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

      typedef typename DiscreteFunctionType :: ConstDofIteratorType
        ConstDofIteratorType;

      //! constructor taking space
      ParallelScalarProduct ( const DiscreteFunctionSpaceType& )
      {}

      ParallelScalarProduct( const ThisType& ) = delete;

      //! return scalar product of dofs
      template < class OtherDiscreteFunctionType >
      RangeFieldType scalarProductDofs ( const DiscreteFunctionType &x,
                                         const OtherDiscreteFunctionType &y ) const
      {
        RangeFieldType scp = 0;

        auto endit = x.dend ();
        auto xit = x.dbegin ();
        auto yit = y.dbegin();
        for( ; xit != endit; ++xit, ++yit )
          scp += (*xit) * (*yit);
        return scp;
      }
    };
#endif

#if HAVE_DUNE_ISTL
    //! Proxy class to evaluate ScalarProduct
    //! holding SlaveDofs which is singleton per space and mapper
    template< class DiscreteFunctionSpaceImp >
    class ParallelScalarProduct
      < ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp > >
    : public ScalarProduct
      < typename ISTLBlockVectorDiscreteFunction< DiscreteFunctionSpaceImp >
          :: DofStorageType >
    {
      //! discrete function type
      typedef ISTLBlockVectorDiscreteFunction<DiscreteFunctionSpaceImp> DiscreteFunctionType;
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
      typedef SlaveDofsProvider< DiscreteFunctionSpaceType > SlaveDofsProviderType;
      typedef typename SlaveDofsProviderType :: SlaveDofsType SlaveDofsType;

      //! export types
      typedef BlockVectorType domain_type;
      typedef typename BlockVectorType :: block_type :: field_type field_type;

      //! define the category
      enum { category=SolverCategory::sequential };

    protected:
      const DiscreteFunctionSpaceType & space_;
      SlaveDofsProviderType slaveDofProvider_;

      ParallelScalarProduct ( const ThisType &org );
    public:
      //! constructor taking space
      ParallelScalarProduct ( const DiscreteFunctionSpaceType &space )
      : space_( space ), slaveDofProvider_( space )
      {}

      //! remove object comm
      ~ParallelScalarProduct () = default;

      SlaveDofsType& slaveDofs ()
      {
        return slaveDofProvider_.slaveDofs();
      }

      const SlaveDofsType& slaveDofs () const
      {
        return slaveDofProvider_.slaveDofs();
      }

      /*! \brief Dot product of two discrete functions.
        It is assumed that the vectors are consistent on the interior+border
        partition.
       */
      RangeFieldType scalarProductDofs(const DiscreteFunctionType& x,
                                       const DiscreteFunctionType& y) const
      {
        std::cout << "ISTL prod" << std::endl;
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
      virtual typename Dune::FieldTraits< RangeFieldType >::real_type
        norm (const BlockVectorType& x)
      {
        return std::abs( std::sqrt( const_cast<ThisType&> (*this).scalarProductDofs(x,x) ) );
      }

      //! delete slave values (for debugging)
      void deleteNonInterior(BlockVectorType& x) const
      {
#if HAVE_MPI
        // case of ALUGrid and DGSpace or FVSpace
        const bool deleteGhostEntries = (space_.gridPart().grid().overlapSize( 0 ) == 0) && !space_.continuous();

        // only delete ghost entries
        if( deleteGhostEntries )
        {
          // rebuild slave dofs if grid was changed
          SlaveDofsType &slaveDofs = slaveDofProvider_.slaveDofs();

          // don't delete the last since this is the overall Size
          const int slaveSize = slaveDofs.size() - 1;
          for(int slave = 0; slave<slaveSize; ++slave)
            x[ slaveDofs[slave] ] = 0;
        }
#endif
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
        SlaveDofsType &slaveDofs = slaveDofProvider_.slaveDofs();

        RangeFieldType scp = 0;
        int i = 0;
        const int slaveSize = slaveDofs.size();
        for(int slave = 0; slave<slaveSize; ++slave)
        {
          const int nextSlave = slaveDofs[slave];
          for(; i<nextSlave; ++i)
          {
            scp += x[i] * y[i];
          }
          // set i to next valid value
          ++i;
        }
        scp = space_.gridPart().comm().sum( scp );
        return scp;
#else
        // return build-in scalar product
        RangeFieldType scp = x * y;
        scp = space_.gridPart().comm().sum( scp );
        return scp;
#endif
      }
    };
#endif

  //@}

  } // end namespace Fem

} // end namespace Dune
#endif // #ifndef DUNE_FEM_SCALARPRODURCTS_HH
