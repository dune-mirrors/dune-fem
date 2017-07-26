#ifndef DUNE_FEM_SPACE_REMAPPED_SPACE_HH
#define DUNE_FEM_SPACE_REMAPPED_SPACE_HH

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

      typedef SlaveDofs< ThisType, BlockMapperType > SlaveDofsType;
    protected:
      typedef typename SlaveDofsType::SingletonKey SlaveDofsKeyType;
      typedef SingletonList< SlaveDofsKeyType, SlaveDofsType > SlaveDofsProviderType;

      struct SlaveDofsDeleter
      {
        void operator() ( SlaveDofsType *slaveDofs )
        {
          SlaveDofsProviderType::removeObject( *slaveDofs );
        }
      };

    public:
      static const int localBlockSize = Traits::localBlockSize;
      typedef Hybrid::IndexRange< int, localBlockSize > LocalBlockIndices;

      ReMappedDiscreteFunctionSpace ( GridPartType &gridPart,
                                  const InterfaceType commInterface = InteriorBorder_All_Interface,
                                  const CommunicationDirection commDirection = ForwardCommunication )
        : BaseType( gridPart, commInterface, commDirection ), blockMapper_( BaseType::blockMapper() )
      {}

      BlockMapperType &blockMapper () const
      {
        return blockMapper_;
      }

      const SlaveDofsType &slaveDofs () const
      {
        if( !slaveDofs_ )
          slaveDofs_.reset( &( SlaveDofsProviderType :: getObject( SlaveDofsKeyType( BaseType::gridPart(), blockMapper() ) ) ),
               SlaveDofsDeleter() );
        slaveDofs_->rebuild( *this );
        return *slaveDofs_;
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
      mutable BlockMapperType blockMapper_;
      mutable std::unique_ptr< CommunicationManagerType > comm_;
      mutable std::shared_ptr< SlaveDofsType > slaveDofs_;
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
