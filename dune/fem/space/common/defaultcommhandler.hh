#ifndef DUNE_FEM_DEFAULTDATAHANDLE_HH
#define DUNE_FEM_DEFAULTDATAHANDLE_HH

#include <cassert>

//- Dune includes
#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/common/hybrid.hh>
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
        typename DiscreteFunction::DofType >
    {
      typedef DefaultCommunicationHandler< DiscreteFunction, Operation > ThisType;
      typedef CommDataHandleIF< ThisType, typename DiscreteFunction::DofType > BaseType;

    public:
      typedef typename BaseType::DataType DataType;

      typedef DiscreteFunction DiscreteFunctionType;

      typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType DiscreteFunctionSpaceType;

    protected:
      typedef typename DiscreteFunctionSpaceType::BlockMapperType BlockMapperType;
      typedef typename DiscreteFunctionSpaceType::LocalBlockIndices LocalBlockIndices;

    public:
      DefaultCommunicationHandler( DiscreteFunctionType &function, const Operation& operation = Operation() )
      : function_( &function ),
        blockMapper_( function.space().blockMapper() ),
        operation_( operation )
      {}

      DefaultCommunicationHandler( const DefaultCommunicationHandler &other )
      : function_( other.function_ ),
        blockMapper_( other.blockMapper_ ),
        operation_( other.operation_ )
      {}

      //! cannot be implemented because of the reference
      DefaultCommunicationHandler &operator= ( const DefaultCommunicationHandler & ) = delete;

    private:
      template < class Buffer >
      struct GatherFunctor
      {
        Buffer& buffer_;
        DiscreteFunctionType *const function_;

        GatherFunctor( Buffer& buffer, DiscreteFunctionType* function )
          : buffer_( buffer ), function_( function )
        {}

        template <class GlobalKey>
        void operator () ( const size_t local, const GlobalKey& globalKey ) const
        {
          const auto &block = function_->dofVector()[ globalKey ];
          Hybrid::forEach( LocalBlockIndices(), [ this, &block ] ( auto &&j ) { buffer_.write( block[ j ] ); } );
        }
      };

      template < class Buffer >
      struct ScatterFunctor
      {
        Buffer& buffer_;
        DiscreteFunctionType *const function_;
        const Operation& operation_;

        ScatterFunctor( Buffer& buffer, DiscreteFunctionType* function, const Operation& operation )
          : buffer_( buffer ), function_( function ), operation_( operation )
        {}

        template <class GlobalKey>
        void operator () ( const size_t local, const GlobalKey& globalKey ) const
        {
          auto &&block = function_->dofVector()[ globalKey ];
          Hybrid::forEach( LocalBlockIndices(), [ this, &block ] ( auto &&j ) {
              DataType value;
              buffer_.read( value );
              operation_( value, block[ j ] );
            } );
        }
      };

    public:
      bool contains ( int dim, int codim ) const
      {
        return blockMapper_.contains( codim );
      }

      bool fixedSize ( int dim, int codim) const
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
        assert( n == size( entity ) );
        ScatterFunctor< MessageBuffer > scatterDofs ( buffer, function_, operation_ );

        blockMapper_.mapEachEntityDof( entity, scatterDofs );
      }

      //! return local dof size to be communicated
      template< class Entity >
      size_t size ( const Entity &entity ) const
      {
        return Hybrid::size( LocalBlockIndices() ) * blockMapper_.numEntityDofs( entity );
      }

    protected:
      DiscreteFunctionType *const function_;
      const BlockMapperType &blockMapper_;
      const Operation operation_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_DEFAULTDATAHANDLE_HH
