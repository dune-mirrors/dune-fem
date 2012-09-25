#ifndef DUNE_FEM_DEFAULTDATAHANDLE_HH
#define DUNE_FEM_DEFAULTDATAHANDLE_HH

#include <cassert>

//- Dune includes 
#include <dune/grid/common/datahandleif.hh>
#include <dune/fem/space/common/commoperations.hh>

namespace Dune
{

  namespace Fem 
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
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;

      typedef typename DiscreteFunctionType::DofBlockPtrType DofBlockPtrType;

      static const int blockSize = DiscreteFunctionSpaceType::localBlockSize;
      
    public:
      DefaultCommunicationHandler( DiscreteFunctionType &function )
      : function_( &function ),
        blockMapper_( function.space().blockMapper() )
      {} 
      
      DefaultCommunicationHandler( const DefaultCommunicationHandler &other )
      : function_( other.function_ ),
        blockMapper_( other.blockMapper_ )
      {}
      
    private:  
      //! cannot be implemented because of the reference
      DefaultCommunicationHandler &operator= ( const DefaultCommunicationHandler & );

      template < class Buffer >
      struct GatherFunctor
      {
        Buffer& buffer_;
        DiscreteFunctionType *const function_;

        GatherFunctor( Buffer& buffer, DiscreteFunctionType* function )
          : buffer_( buffer ),
            function_( function )
        {
        }

        template <class GlobalKey>
        void operator () ( const size_t local, const GlobalKey& globalKey ) 
        {
          DofBlockPtrType blockPtr = function_->block( globalKey );
          for( int j = 0; j < blockSize; ++j )
          {
            buffer_.write( (*blockPtr)[ j ] );
          }
        }
      };

      template < class Buffer >
      struct ScatterFunctor
      {
        Buffer& buffer_;
        DiscreteFunctionType *const function_;

        ScatterFunctor( Buffer& buffer, DiscreteFunctionType* function )
          : buffer_( buffer ),
            function_( function )
        {
        }

        template <class GlobalKey>
        void operator () ( const size_t local, const GlobalKey& globalKey ) 
        {
          DofBlockPtrType blockPtr = function_->block( globalKey );
          for( int j = 0; j < blockSize; ++j )
          {
            DataType value;
            buffer_.read( value );

            Operation :: apply( value, (*blockPtr)[ j ] );
          }
        }
      };
    public:
      bool contains ( int dim, int codim ) const
      {
        return blockMapper_.contains( codim );
      }

      bool fixedsize ( int dim, int codim) const
      {
        return blockMapper_.fixedDataSize( codim );
      }

      //! read buffer and apply operation 
      template< class MessageBuffer, class Entity >
      void gather ( MessageBuffer &buffer, const Entity &entity ) const
      {
        GatherFunctor< MessageBuffer > gatherDofs ( buffer, function_ );
        blockMapper_.mapEachEntityDof( entity, gatherDofs );
      }

      //! read buffer and apply operation 
      template< class MessageBuffer, class Entity >
      void scatter ( MessageBuffer &buffer, const Entity &entity, size_t n )
      {
        assert( n == blockSize *  blockMapper_.numEntityDofs( entity ) );
        ScatterFunctor< MessageBuffer > scatterDofs ( buffer, function_ );

        blockMapper_.mapEachEntityDof( entity, scatterDofs );
      }

      //! return local dof size to be communicated 
      template< class Entity >
      size_t size ( const Entity &entity ) const
      {
        return blockSize * blockMapper_.numEntityDofs( entity );
      }

    protected:
      DiscreteFunctionType *const function_;
      const BlockMapperType &blockMapper_;
    };
  
  } // namespace Fem 

} // namespace Dune

#endif // #ifndef DUNE_FEM_DEFAULTDATAHANDLE_HH
