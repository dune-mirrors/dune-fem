#ifndef DUNE_FEM_COMMOPERATIONS_HH
#define DUNE_FEM_COMMOPERATIONS_HH

#include <tuple>
#include <type_traits>
#include <utility>

#include <dune/common/tupleutility.hh>

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

      /** DataGather functor  */
      template <class BufferImp, class EntityImp>
      class DataGather{
      public:
        //! Constructor taking buffer and entity
        DataGather(BufferImp & buff, const EntityImp & en )
        : buff_(buff)
        , en_(en)
        {}

        //! call gather on given data handle object
        template <class DataHandlerImp>
        void visit(DataHandlerImp & dh) {
          dh.gather(buff_,en_);
        }

      private:
        BufferImp & buff_;
        const EntityImp & en_;
      };

      /** DataScatter functor  */
      template <class BufferImp, class EntityImp>
      class DataScatter{
      public:
        //! Constructor
        //! Constructor taking buffer and entity and size
        DataScatter(BufferImp & buff, const EntityImp & en, size_t n)
        : buff_(buff)
        , en_(en)
        , size_(n)
        {}

        //! call scatter on given data handle object
        template <class DataHandlerImp>
        void visit(DataHandlerImp & dh) {
          // TODO: here, the wrong size is passed to the subhandles
          dh.scatter(buff_,en_,size_);
        }

      private:
        BufferImp & buff_;
        const EntityImp & en_;
        const size_t size_;
      };

      /** DataSize functor  */
      template <class EntityImp>
      class DataSize
      {
      public:
        //! Constructor
        //! \param entity to calc size for
        DataSize(const EntityImp & en)
        : en_(en)
        , size_(0)
        {}

        //! call size on given data handle object
        template <class DataHandlerImp>
        void visit(DataHandlerImp & dh)
        {
          size_ += dh.size(en_);
        }

        //! return size
        size_t size() const { return size_; }

      private:
        const EntityImp & en_;
        size_t size_;
      };

      /** FixedSize functor  */
      class FixedSize
      {
      public:
        //! Constructor
        //! \param dim to check for
        //! \param codim to check for
        FixedSize(const int dim, const int codim)
        : dim_(dim)
        , codim_(codim)
        , fixedSize_(true)
        {}

        //! call size on given data handle object
        template <class DataHandlerImp>
        void visit(DataHandlerImp & dh)
        {
          bool fs = dh.fixedsize(dim_,codim_);
          fixedSize_ = (fs == false) ? fs : fixedSize_;
        }

        //! return size
        bool fixedSize() const { return fixedSize_; }

      private:
        const int dim_;
        const int codim_;
        bool fixedSize_;
      };

      /** Contains functor  */
      class Contains
      {
      public:
        //! Constructor
        //! \param dim to check for
        //! \param codim to check for
        Contains(const int dim, const int codim)
        : dim_(dim)
        , codim_(codim)
        , contains_(false)
        {}

        //! call size on given data handle object
        template <class DataHandlerImp>
        void visit(DataHandlerImp & dh)
        {
          bool c = dh.contains(dim_,codim_);
          contains_ = (c == true) ? c : contains_;
        }

        //! return size
        bool contains() const { return contains_; }

      private:
        const int dim_;
        const int codim_;
        bool contains_;
      };

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
        ForEachValue<DataHandlerTupleType> forEach(data_);
        Contains dataContains(dim,codim);
        forEach.apply(dataContains);
        return dataContains.contains();
      }

      bool fixedsize (int dim, int codim) const
      {
        ForEachValue<DataHandlerTupleType> forEach(data_);
        FixedSize dataFixedSize(dim,codim);
        forEach.apply(dataFixedSize);
        return dataFixedSize.fixedSize();
      }

      //! \brief loop over all internal data handlers and call gather for
      //! given entity
      template<class MessageBufferImp, class EntityType>
      void gather (MessageBufferImp& buff, const EntityType& en) const
      {
        ForEachValue<DataHandlerTupleType> forEach(data_);
        DataGather<MessageBufferImp,EntityType> gatherData(buff,en);
        forEach.apply(gatherData);
      }

      //! \brief loop over all internal data handlers and call scatter for
      //! given entity
      template<class MessageBufferImp, class EntityType>
      void scatter (MessageBufferImp& buff, const EntityType& en, size_t n)
      {
        ForEachValue<DataHandlerTupleType> forEach(data_);
        DataScatter<MessageBufferImp,EntityType> scatterData(buff,en,n);
        forEach.apply(scatterData);
      }

      //! \brief loop over all internal data handlers and return sum of data
      //! size of given entity
      template<class EntityType>
      size_t size (const EntityType& en) const
      {
        ForEachValue<DataHandlerTupleType> forEach(data_);
        DataSize<EntityType> dataSize(en);
        forEach.apply(dataSize);
        return dataSize.size();
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
      enum dfCommunicationOperation {copy,add,sub,min,max};
      //! just copy data
      struct Copy
      {
        static const dfCommunicationOperation value = copy;
        static const char * name ()
        {
          return "Copy";
        }

        template <class DataType>
        static inline void apply(const DataType & arg, DataType & dest)
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
        static inline void apply(const DataType & arg, DataType & dest)
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
        static inline void apply(const DataType & arg, DataType & dest)
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
        static inline void apply(const DataType & arg, DataType & dest)
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
        static inline void apply(const DataType & arg, DataType & dest)
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
