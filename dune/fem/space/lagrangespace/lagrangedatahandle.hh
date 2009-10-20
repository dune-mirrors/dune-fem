#ifndef DUNE_LAGRANGESPACEDATAHANDLE_HH
#define DUNE_LAGRANGESPACEDATAHANDLE_HH

#include <cassert>

//- Dune includes 
#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/space/common/commoperations.hh>

namespace Dune
{

  /** \class   LagrangeCommunicationHandler
   *  \ingroup DFComm
   *  \brief   Communication handler for Lagrange discrete functions
   *
   *  \param  DiscreteFunction  type of discrete function to be communicated
   */
  template< class DiscreteFunction,
            class Operation = DFCommunicationOperation :: Add >
  class LagrangeCommunicationHandler
  : public CommDataHandleIF
    < LagrangeCommunicationHandler< DiscreteFunction, Operation >,
      typename DiscreteFunction :: RangeFieldType >
  {
  public:  
    typedef DiscreteFunction DiscreteFunctionType;
    
    typedef typename DiscreteFunctionType :: RangeFieldType DataType;

  private:
    typedef CommDataHandleIF< LagrangeCommunicationHandler, DataType > BaseType;

  public:
    typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType
      DiscreteFunctionSpaceType;

    typedef typename DiscreteFunctionSpaceType :: BlockMapperType MapperType;

  protected:
    typedef typename DiscreteFunctionType :: DofBlockPtrType DofBlockPtrType;

    enum { blockSize = DiscreteFunctionSpaceType :: localBlockSize };
    
  protected:
    DiscreteFunctionType *const function_;
    const MapperType &mapper_;

  public:
    LagrangeCommunicationHandler( DiscreteFunctionType &function )
    : function_( &function ),
      mapper_( function.space().blockMapper() )
    {} 
    
    LagrangeCommunicationHandler( const LagrangeCommunicationHandler &other )
    : function_( other.function_ ),
      mapper_( other.mapper_ )
    {}
    
  private:  
    //! cannot be implemented because of the reference
    LagrangeCommunicationHandler &operator= ( const LagrangeCommunicationHandler & );

  public:
    inline bool contains ( int dim, int codim ) const
    {
      return mapper_.contains( codim );
    }

    inline bool fixedsize ( int dim, int codim) const
    {
      return mapper_.fixedDataSize( codim );
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void gather ( MessageBuffer &buffer, const Entity &entity ) const
    {
      const unsigned int numEntityDofs = mapper_.numEntityDofs( entity );
      for( unsigned int i = 0; i < numEntityDofs; ++i )
      {
        const unsigned int index = mapper_.mapEntityDofToGlobal( entity, i );
        
        DofBlockPtrType blockPtr = function_->block( index );
        for( unsigned int j = 0; j < blockSize; ++j )
          buffer.write( (*blockPtr)[ j ] );
      }
    }

    //! read buffer and apply operation 
    template< class MessageBuffer, class Entity >
    void scatter ( MessageBuffer &buffer, const Entity &entity, size_t n )
    {
      const unsigned int numEntityDofs = mapper_.numEntityDofs( entity );

      assert( n == numEntityDofs );
      for( unsigned int i = 0; i < numEntityDofs; ++i )
      {
        const unsigned int index = mapper_.mapEntityDofToGlobal( entity, i );

        DofBlockPtrType blockPtr = function_->block( index );
        for( unsigned int j = 0; j < blockSize; ++j )
        {
          DataType value;
          buffer.read( value );

          Operation :: apply( value, (*blockPtr)[ j ] );
        }
      }
    }

    //! return local dof size to be communicated 
    template< class Entity >
    size_t size ( const Entity &entity ) const
    {
      return blockSize * mapper_.numEntityDofs( entity );
    }
  };
  
} // end namespace Dune
#endif
