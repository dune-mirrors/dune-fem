#ifndef DUNE_FEM_COMMOPERATIONS_HH
#define DUNE_FEM_COMMOPERATIONS_HH

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/grid/common/datahandleif.hh>

namespace Dune
{

  namespace Fem
  {

  /** @addtogroup DFComm
      @{
  **/

    // CombinedDataType
    // ----------------

    template< class... DataHandle >
    struct CombinedDataType;

    template< class DataHandle >
    struct CombinedDataType< DataHandle >
    {
      typedef typename DataHandle::DataType Type;
    };

    template< class DataHandle, class... Tail >
    struct CombinedDataType< DataHandle, Tail... >
    {
      typedef typename DataHandle::DataType Type;
      static_assert( std::is_same< Type, typename CombinedDataType< Tail... >::Type >::value, "Only data handles for the same data type can be combined." );
    };



    // CombinedDataHandle
    // ------------------

    /**
     * \brief combine multiple data handles into one
     */
    template< class... DataHandle >
    class CombinedDataHandle
      : public CommDataHandleIF< CombinedDataHandle< DataHandle... >, typename CombinedDataType< DataHandle... >::Type >
    {
      typedef std::tuple< DataHandle... > DataHandlerTupleType;
      static constexpr std::size_t tupleSize = std::tuple_size< DataHandlerTupleType >::value;

    public:
      typedef typename CombinedDataType< DataHandle... >::Type DataType;

      CombinedDataHandle ( const DataHandle &... handle )
        : data_( handle... )
      {}

      CombinedDataHandle ( const std::tuple< DataHandle... > &data )
        : data_( data )
      {}

      bool contains (int dim, int codim) const
      {
        bool value( false );
        Hybrid::forEach( std::make_index_sequence< tupleSize >{},
          [ & ]( auto i ){ value = ( value || std::get< i >( data_ ).contains( dim, codim ) ); } );
        return value;
      }

      bool fixedSize (int dim, int codim) const
      {
        bool value( true );
        Hybrid::forEach( std::make_index_sequence< tupleSize >{},
          [ & ]( auto i ){ value = ( value && std::get< i >( data_ ).fixedSize( dim, codim ) ); } );
        return value;
      }

      //! \brief loop over all internal data handlers and call gather for
      //! given entity
      template<class MessageBufferImp, class EntityType>
      void gather (MessageBufferImp& buff, const EntityType& en) const
      {
        Hybrid::forEach( std::make_index_sequence< tupleSize >{}, [ & ]( auto i ){ std::get< i >( data_ ).gather( buff, en ); } );
      }

      //! \brief loop over all internal data handlers and call scatter for
      //! given entity
      template<class MessageBufferImp, class EntityType>
      void scatter (MessageBufferImp& buff, const EntityType& en, std::size_t n)
      {
        Hybrid::forEach( std::make_index_sequence< tupleSize >{}, [ & ]( auto i ){ std::get< i >( data_ ).scatter( buff, en, n ); } );
      }

      //! \brief loop over all internal data handlers and return sum of data
      //! size of given entity
      template<class EntityType>
      std::size_t size (const EntityType& en) const
      {
        std::size_t value( 0 );
        Hybrid::forEach( std::make_index_sequence< tupleSize >{}, [ & ]( auto i ){ value += std::get< i >( data_ ).size( en ); } );
        return value;
      }

    private:
      DataHandlerTupleType data_;
    };

    ////////////////////////////////////////////////////////////////
    //
    //  --DiscreteFunctionCommunications
    //
    ////////////////////////////////////////////////////////////////

    //! \brief Mathematical operation apply during communication
    //! to data that is communicated
    //! enum of all avialable operations
    struct DFCommunicationOperation
    {
      enum dfCommunicationOperation { copy, add, sub, min, max };

      //! just copy data
      struct Copy
      {
        static const dfCommunicationOperation value = copy;
        static const char * name ()
        {
          return "Copy";
        }

        template <class DataType>
        inline void operator () (const DataType & arg, DataType & dest) const
        {
          dest = arg;
        }
      };


      //! sum up data
      struct Add
      {
        static const dfCommunicationOperation value = add;
        static const char * name ()
        {
          return "Add";
        }

        template <class DataType>
        inline void operator () (const DataType & arg, DataType & dest) const
        {
          dest += arg;
        }
      };

      //! substract data
      struct Sub
      {
        static const dfCommunicationOperation value = sub;
        static const char * name ()
        {
          return "Sub";
        }

        template <class DataType>
        inline void operator () (const DataType & arg, DataType & dest) const
        {
          dest -= arg;
        }
      };

      //! keep minimum
      struct Min
      {
        static const dfCommunicationOperation value = min;
        static const char * name ()
        {
          return "Min";
        }

        template <class DataType>
        inline void operator () (const DataType & arg, DataType & dest) const
        {
          dest = std::min(dest,arg);
        }
      };

      //! keep maximum
      struct Max
      {
        static const dfCommunicationOperation value = max;
        static const char * name ()
        {
          return "Max";
        }

        template <class DataType>
        inline void operator () (const DataType & arg, DataType & dest) const
        {
          dest = std::max(dest,arg);
        }
      };
    };

    ////////////////////////////////////////////////////////////////
    //
    //  --LoadBalanceContainsCheck
    //
    ////////////////////////////////////////////////////////////////

    //! \brief check for sets of entities for the load balance procedure
    template <class DiscreteFunction>
    class LoadBalanceLeafData
    {
    public:
      typedef DiscreteFunction DiscreteFunctionType ;
      typedef typename DiscreteFunctionType :: DiscreteFunctionSpaceType :: IteratorType :: Entity Entity;

      explicit LoadBalanceLeafData( const DiscreteFunctionType& df ) {}
      /** \brief return true if the data of this entity should be transfered during load balance */
      bool contains (const Entity& entity) const
      {
        return entity.isLeaf();
      }
    };

  //@}

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMOPERATIONS_HH
