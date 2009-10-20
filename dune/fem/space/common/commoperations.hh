#ifndef DUNE_COMMOPERATIONS_HH
#define DUNE_COMMOPERATIONS_HH

#include <dune/fem/misc/utility.hh>
//- Dune includes 

#include <dune/grid/common/datahandleif.hh>

namespace Dune {

/** @addtogroup DFComm  
    @{
**/

  /** Class to combine data handles. The outcome is a class satifying the
      DataHandleIF itself. 
   */
  template <  class DataImpOne 
            , class DataImpTwo   = Nil  
            , class DataImpThree = Nil
            , class DataImpFour  = Nil 
            , class DataImpFive  = Nil 
            , class DataImpSix   = Nil  
            , class DataImpSeven = Nil  
            , class DataImpEight = Nil  
            , class DataImpNine  = Nil > 
  class CombinedDataHandle 
   : public CommDataHandleIF< 
       CombinedDataHandle < DataImpOne , 
                            DataImpTwo ,
                            DataImpThree,
                            DataImpFour ,
                            DataImpFive ,
                            DataImpSix  , 
                            DataImpSeven, 
                            DataImpEight, 
                            DataImpNine 
                            > , 
                            typename DataImpOne ::DataType >
  {
    /** DataGather functor  */
    template <class BufferImp, class EntityImp>
    class DataGather{
    public:
      //! Constructor taking buffer and entity 
      DataGather(BufferImp & buff, const EntityImp & en) 
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

    DataImpOne   & one_;
    DataImpTwo   & two_;
    DataImpThree & three_;
    DataImpFour  & four_;
    DataImpFive  & five_;
    DataImpSix   & six_;
    DataImpSeven & seven_;
    DataImpEight & eight_;
    DataImpNine  & nine_;

    typedef Tuple<  DataImpOne 
                  , DataImpTwo 
                  , DataImpThree 
                  , DataImpFour 
                  , DataImpFive  
                  , DataImpSix
                  , DataImpSeven 
                  , DataImpEight
                  , DataImpNine 
                  > DataHandlerTupleType;

    mutable DataHandlerTupleType data_; 
  public:
    typedef typename DataImpOne :: DataType DataType;

    // initialize of Nil 
    static Nil & null() { 
      static Nil n;
      return n;
    }

    CombinedDataHandle(DataImpOne & one
              , DataImpTwo   & two   = null()  
              , DataImpThree & three = null()
              , DataImpFour  & four  = null()
              , DataImpFive  & five  = null() 
              , DataImpSix   & six   = null() 
              , DataImpSeven & seven = null() 
              , DataImpEight & eight = null() 
              , DataImpNine  & nine  = null() )
      : one_(one) , two_(two) , three_(three) , four_(four) 
      , five_(five) , six_(six_) , seven_(seven)
      , eight_(eight) , nine_(nine_) 
      , data_( one_ , two_, three_, four_, five_, six_ , seven_ , eight_ , nine_ )
    {
    }

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
  };

  ////////////////////////////////////////////////////////////////
  //
  //  --DiscreteFunctionCommunications 
  //
  ////////////////////////////////////////////////////////////////

  //! \brief Mathematical operation apply during communication 
  //! to data that is communicated 
  struct DFCommunicationOperation 
  {
    //! just copy data 
    struct Copy 
    {
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
    
    //! keep minimum   
    struct Min 
    {
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
 
//@} 
} // end namespace Dune
#endif
