#ifndef DUNE_FEM_DEFAULTDATAHANDLE_HH
#define DUNE_FEM_DEFAULTDATAHANDLE_HH

#include <cassert>

//- Dune includes 
#include <dune/grid/common/datahandleif.hh>
#include <dune/fem/space/common/commoperations.hh>

namespace Dune
{

  /** \class   DefaultCommunicationHandler
   *  \ingroup DFComm
   *  \brief   Default communication handler for discrete functions
   *
   *  \param  DiscreteFunction  type of discrete function to be communicated
   */
  template< class DiscreteFunction, class Operation = DFCommunicationOperation::Add >
  class DefaultCommunicationHandler
  : public CommDataHandleIF
    < DefaultCommunicationHandler< DiscreteFunction, Operation >,
      typename DiscreteFunction::RangeFieldType >
  {
    typedef DefaultCommunicationHandler< DiscreteFunction, Operation > ThisType;
    typedef CommDataHandleIF< ThisType, typename DiscreteFunction::RangeFieldType > BaseType;

  public:  
    typedef typename BaseType::DataType DataType;

    typedef DiscreteFunction DiscreteFunctionType;

    typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

  protected:
    typedef typename DiscreteFunctionSpaceType::BlockMapperType MapperType;

    typedef typename DiscreteFunctionType::DofBlockPtrType DofBlockPtrType;

    static const int blockSize = DiscreteFunctionSpaceType::localBlockSize;
    
  public:
    DefaultCommunicationHandler( DiscreteFunctionType &function )
    : function_( &function ),
      mapper_( function.space().blockMapper() )
    {} 
    
    DefaultCommunicationHandler( const DefaultCommunicationHandler &other )
    : function_( other.function_ ),
      mapper_( other.mapper_ )
    {}
    
  private:  
    //! cannot be implemented because of the reference
    DefaultCommunicationHandler &operator= ( const DefaultCommunicationHandler & );

  public:
    bool contains ( int dim, int codim ) const
    {
      return mapper_.contains( codim );
    }

    bool fixedsize ( int dim, int codim) const
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
        for( int j = 0; j < blockSize; ++j )
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
        for( int j = 0; j < blockSize; ++j )
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

  protected:
    DiscreteFunctionType *const function_;
    const MapperType &mapper_;
  };
  
} // end namespace Dune

#endif // #ifndef DUNE_FEM_DEFAULTDATAHANDLE_HH
