#ifndef DUNE_FEM_SPACE_REMAPPED_SPACE_HH
#define DUNE_FEM_SPACE_REMAPPED_SPACE_HH

#include <memory>
#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/fem/common/hybrid.hh>
#include <dune/fem/space/common/discretefunctionspace.hh>

namespace Dune
{

  namespace Fem
  {

    // internal forward declaration

    // ReMappedDiscreteFunctionSpace
    // -----------------------------

    template< class DiscreteFunctionSpace, class BlockMapper, int blockSize >
    class ReMappedDiscreteFunctionSpace;


    // ReMappedDiscreteFunctionSpaceTraits
    // -----------------------------------

    template< class DiscreteFunctionSpace, class BlockMapper, int blockSize >
    struct ReMappedDiscreteFunctionSpaceTraits
      : public DiscreteFunctionSpace::Traits
    {
      typedef ReMappedDiscreteFunctionSpace< DiscreteFunctionSpace, BlockMapper, blockSize > DiscreteFunctionSpaceType;

      static const int localBlockSize = blockSize;

      typedef BlockMapper BlockMapperType;
    };


    // ReMappedDiscreteFunctionSpace
    // -----------------------------

    template< class DiscreteFunctionSpace, class BlockMapper, int blockSize >
    class ReMappedDiscreteFunctionSpace
      : public DiscreteFunctionSpace
    {
      typedef ReMappedDiscreteFunctionSpace< DiscreteFunctionSpace, BlockMapper, blockSize > ThisType;
      typedef DiscreteFunctionSpace BaseType;

    public:
      typedef ReMappedDiscreteFunctionSpaceTraits< DiscreteFunctionSpace, BlockMapper, blockSize > Traits;

    private:
      // here we need to declare the BaseType::{blockMapper(), communicator(), slaveDofs()} private
      using BaseType::blockMapper;
      using BaseType::communicator;
      using BaseType::slaveDofs;

      typedef CommunicationManager< ThisType > CommunicationManagerType;
    public:

      typedef typename BaseType::GridPartType GridPartType;
      typedef typename Traits::BlockMapperType BlockMapperType;

      template< class DiscreteFunction, class Operation >
      using CommDataHandleType = typename BaseType::template CommDataHandle< DiscreteFunction, Operation >::Type;

      typedef SlaveDofs< ThisType, BlockMapperType > SlaveDofsType;
    protected:
      struct SlaveDofsFactory
      {
        typedef std::pair< SlaveDofsType, int > ObjectType;

        static ObjectType *createObject ( std::pair< GridPartType *, BlockMapperType * > key )
        {
          return new ObjectType( std::piecewise_construct, std::tie( *key.first, *key.second ), std::make_tuple( -1 ) );
        }

        static void deleteObject ( ObjectType *object ) { delete object; }
      };

      typedef SingletonList< std::pair< GridPartType *, BlockMapperType * >, std::pair< SlaveDofsType, int >, SlaveDofsFactory > SlaveDofsProviderType;

    public:
      static const int localBlockSize = Traits::localBlockSize;
      typedef Hybrid::IndexRange< int, localBlockSize > LocalBlockIndices;

      ReMappedDiscreteFunctionSpace ( GridPartType &gridPart,
                                  const InterfaceType commInterface = InteriorBorder_All_Interface,
                                  const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection )
      {
        createBlockMapper( gridPart, BaseType::BlockMapper(), blockMapper_ );
      }

      BlockMapperType &blockMapper () const
      {
        assert( blockMapper_ );
        return blockMapper_;
      }

      const SlaveDofsType &slaveDofs () const
      {
        if( !slaveDofs_ )
          slaveDofs_.reset( &( SlaveDofsProviderType :: getObject( std::make_pair( &BaseType::gridPart(), &blockMapper() ) ) ) );
        const int sequence = this->sequence();
        if( slaveDofs_->second != sequence )
        {
          slaveDofs_->first.rebuild();
          slaveDofs_->second = sequence;
        }
        return slaveDofs_->first;
      }

      template< class DiscreteFunction, class Operation >
      CommDataHandleType< DiscreteFunction, Operation >
      createDataHandle( DiscreteFunction &discreteFunction, const Operation &operation ) const
      {
        static_assert( std::is_same< BlockMapperType, typename DiscreteFunction::DiscreteFunctionSpaceType::BlockMapperType >::value
                       && (localBlockSize == static_cast< std::size_t >( DiscreteFunction::DiscreteFunctionSpaceType::localBlockSize )),
                       "ReMappedDiscreteFunctionSpace::createDataHandle cannot be called with discrete functions defined over a different space" );
        return CommDataHandleType< DiscreteFunction, Operation >( discreteFunction, operation );
      }

      const CommunicationManagerType &communicator () const
      {
        if( !comm_ )
          comm_.reset( new CommunicationManagerType( *this, BaseType::communicationInterface(), BaseType::communicationDirection() ) );
        return *comm_;
      }

      template< class DiscreteFunction >
      void communicate ( DiscreteFunction &df ) const
      {
        typename Traits::template CommDataHandle< DiscreteFunction >::OperationType op;
        communicate( df, op );
      }

      template< class DiscreteFunction, class Op >
      void communicate ( DiscreteFunction &df, const Op &op ) const
      {
        communicator().exchange( df, op );
      }

    private:
      template< class BM >
      static std::enable_if_t< std::is_constructible< BM, typename BaseType::BlockMapperType & >::value >
      createBlockMapper ( const GridPartType &gridPart, typename BaseType::BlockMapperType &baseMapper, std::unique_ptr< BM > &blockMapper )
      {
        blockMapper.reset( new BM( baseMapper ) );
      }

      template< class BM >
      static std::enable_if_t< std::is_constructible< BM, const GridPartType &, typename BaseType::BlockMapperType & >::value >
      createBlockMapper ( const GridPartType &gridPart, typename BaseType::BlockMapperType &baseMapper, std::unique_ptr< BM > &blockMapper )
      {
        blockMapper.reset( new BM( gridPart, baseMapper ) );
      }

      std::unique_ptr< BlockMapperType > blockMapper_;
      mutable std::unique_ptr< CommunicationManagerType > comm_;
      mutable std::unique_ptr< std::pair< SlaveDofsType, int >, typename SlaveDofsProviderType::Deleter > slaveDofs_;
    };


    // DifferentDiscreteFunctionSpace
    // ------------------------------

    template< class DiscreteFunctionSpace, class BlockMapper, int blockSize, class NewFunctionSpace >
    struct DifferentDiscreteFunctionSpace< ReMappedDiscreteFunctionSpace< DiscreteFunctionSpace, BlockMapper, blockSize >, NewFunctionSpace >
    {
      typedef typename DifferentDiscreteFunctionSpace< DiscreteFunctionSpace, NewFunctionSpace >::Type NewSpaceType;

      static const int newBlockSize = blockSize * NewSpaceType::localBlockSize / DiscreteFunctionSpace::localBlockSize;
      typedef ReMappedDiscreteFunctionSpace< NewSpaceType, BlockMapper, newBlockSize > Type;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_SPACE_REMAPPED_SPACE_HH
