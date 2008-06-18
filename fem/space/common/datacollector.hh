#ifndef DUNE_DATACOLLECTOR_HH
#define DUNE_DATACOLLECTOR_HH

//-System includes
#include <cassert>
#include <vector>
#include <utility>
#include <iostream>

//-Dune includes 
#include <dune/common/interfaces.hh>

//- local includes 
#include <dune/fem/operator/common/objpointer.hh>

namespace Dune{

/** @addtogroup DataCollectors 
    
  @{
 */

/*! This could be seen as a hack, but is not
  With this method we define the class CombineRestProl which is only 
  for combining of local grid operations without using virtual methods.
  The interface of the two defined methods of the class (restrictLocal and
  prolongLocal) is given by the implementation (see below ) and 
  has to be the same for all local operators you want to combine 
*/

// forward declareation of DofManager
template <class GridImp> class DofManager;

template <class LocalOp, class ParamType> class LocalInlinePlus;

#define PARAM_CLASSNAME CombinedLocalDataCollect 
#define PARAM_INHERIT LocalInlinePlus  
#define PARAM_FUNC_1 apply 
#include <dune/fem/operator/common/combine.inc>

template <class ParamT>
class LocalInterface : public ObjPointerStorage
{
public:
  typedef LocalInterface<ParamT> MyType;
  struct Traits 
  {
    typedef ParamT ParamType;
  };

private:
  typedef ParamT ParamType;
  typedef void FuncType(MyType &, ParamType & p); 
  typedef typename std::pair < MyType * ,  FuncType * > PairType;
  typedef typename std::vector < PairType > ListType;
  
  template <class OpType> 
  struct AddToWrapper 
  {
    //! applyWrapper knows the real type of Op 
    static void applyWrapper(MyType & m, ParamType & p )
    {
      static_cast<OpType &> (m).apply(p); 
    }

    //! store pointer to Op /as a pointer to interface class 
    //! and store a pointer to applyWrapper as a function pointer
    static void addToList (ListType & vec , const OpType & op )
    {
      PairType p( const_cast<OpType *> (&op) , applyWrapper);
      vec.push_back(p);
    }
  };

  // copy list of op to this class 
  static void copyList (ListType & vec , const MyType & op )
  {
    const ListType & ve = op.vec_;
    if(ve.size() > 0)
    {   
      ListType tmp ( vec );
      vec.resize( vec.size() + ve.size() );
    
      // copy list to new vector 
      for(unsigned int i=0; i<tmp.size(); i++)
        vec[i] = tmp[i];
      for(unsigned int i=tmp.size(); i<vec.size(); i++)
        vec[i] = ve[i-tmp.size()];
    }
  }
  
public:
  LocalInterface () : vec_ (0) {}
  
  template <class OpType>
  LocalInterface (const OpType & op)
  {
    AddToWrapper<OpType>::addToList(vec_,op);
  }

  LocalInterface (const MyType & op)
  {
    copyList(vec_,op);
  }

  virtual ~LocalInterface() {};

  //! for all pointer to local operators call the func pointer 
  void apply ( ParamType & p ) const 
  {
    const size_t size = vec_.size();
    for(size_t i=0; i<size; ++i)
    {
      assert( vec_[i].second );
      assert( vec_[i].first );
      // vec_[i].second contains the pointer to the function that makes the
      // correct cast to the real type of vec_[i].first 
      (*vec_[i].second)( *(vec_[i].first) , p );
    }
  }
    
  template <class OpType>
  MyType & operator + (const OpType & op)
  {
    AddToWrapper<OpType>::addToList(vec_,op);
    return *this;
  }
  
  MyType & operator + (const MyType & op)
  {
    copyList(vec_,op);
    return *this;
  }
  
  template <class OpType>
  MyType & operator += (const OpType & op)
  {
    AddToWrapper<OpType>::addToList(vec_,op);
    return *this;
  }
  
  MyType & operator += (const MyType & op)
  {
    copyList(vec_,op);
    return *this;
  }

  template <class OpType>
  void remove(const OpType & op)
  {
    typedef typename ListType :: iterator iterator; 
    iterator end = vec_.end(); 
    for(iterator it = vec_.begin(); it != end; ++it)
    {
      if( &op == (*it).first ) 
      {
        vec_.erase( it ); 
        return ;
      }
    }
  }
  
  template <class OpType>
  MyType & operator = (const OpType & op) 
  {
    AddToWrapper<OpType>::addToList(vec_,op);
    return *this;
  }

  bool empty () const { return (vec_.size() == 0); }
 
private:  
  mutable ListType vec_;
};

template <class LocalOp, class ParamT>
class LocalInlinePlus : public LocalInterface<ParamT> 
{
public:
  //! when overload this class ParamType  
  //! and LocalInterfaceType must be redefined 
  struct Traits 
  {
    typedef ParamT ParamType;
    typedef LocalInterface<ParamType> LocalInterfaceType;
  };
  
  template <class B> 
  CombinedLocalDataCollect<LocalOp,B> & operator + (const B & b) 
  {
    std::cout << "operator + of LocalInlinePlus \n";
    typedef CombinedLocalDataCollect<LocalOp,B> CombinedType;
    CombinedType * combo = new CombinedType ( asImp() , b );
    saveObjPointer( combo );   
    return *combo;
  }
  LocalOp & asImp() { return static_cast<LocalOp &> (*this); }
};

class DummyObjectStream
{ 
  public:
    class EOFException {} ;
    template <class T>
    void read (T &) const { assert(false); abort(); }
    template <class T>
    void readObject (T &) { assert(false); abort(); }
    void readObject (int) { assert(false); abort(); }
    void readObject (double) { assert(false); abort(); }
    template <class T>
    void write (const T &) { assert(false);abort(); }
    template <class T>
    void writeObject (T &) { assert(false);abort(); }
    void writeObject (int) { assert(false);abort(); }
    void writeObject (double) { assert(false);abort(); }
};

/*! Combination of different DataCollectors

 This Class is the result of a combination of different
 AdaptationOperators. It is the same principle as with Mapping and
 DiscreteOperatorImp. 
*/ 
template <class GridType, class ObjectStreamImp = DummyObjectStream >
class DataCollectorInterface 
{
  typedef ObjectStreamImp  ObjectStreamType;
  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef DataCollectorInterface<GridType,ObjectStreamImp> MyType;
  typedef std::pair < ObjectStreamType * , 
          const typename GridType:: template Codim<0>::Entity  * >  DataCollectorParamType;
  
public:
  typedef LocalInterface<DataCollectorParamType> LocalInterfaceType;

  struct Traits 
  {
    typedef LocalInterface<DataCollectorParamType> LocalInterfaceType;
  };
  
  //! empty constructor 
  DataCollectorInterface () : dc_ (0) {}
  
  /** \brief Virtual desctructor */
  virtual ~DataCollectorInterface() { clear(); }

  //! all adaptation operators have this method which adapts the
  //! corresponding grid and organizes the restriction prolongation process
  //! of the underlying function spaces
  virtual void apply (ObjectStreamType &str,EntityType & en) const 
  {
    //std::cout << "apply on interface class \n";
    if(dc_) (*dc_).apply(str,en);  
    else 
    {
      std::cerr << "WARNING: apply: did nothing! \n";
    }
  };
  
  virtual const LocalInterfaceType & getLocalInterfaceOp () const 
  {
    if(dc_) 
      return dc_->getLocalInterfaceOp ();
    else 
    {
      std::cerr << "No LocalInterfaceOperator \n";
      assert(false);
      return *(new LocalInterfaceType());
    }
  };
  
  virtual LocalInterfaceType & getLocalInterfaceOp ()
  {
    if(dc_) 
      return dc_->getLocalInterfaceOp ();
    else 
    {
      std::cerr << "No LocalInterfaceOperator \n";
      assert(false);
      return *(new LocalInterfaceType());
    }
  };

  //! Assignement operator
  template <class OpType>
  MyType & operator += (const OpType & dc)
  {
    if(dc_)
    {
      //std::cout << "Operator += with OpType \n";
      dc_ = dcConv_;
      MyType * tmp = const_cast<OpType &> (dc).convert();
      dc_ = &(*dc_).operator += (*tmp);
    }
    else 
    {
      dc_     = const_cast<OpType *> (&dc);
      dcConv_ = const_cast<OpType &> (dc).convert();
    }
    return (*this);
  }
  
  //! Assignement operator
  virtual MyType & operator += (const MyType & dc)
  {
    if(dc_)
    {
      //std::cout << "Operator += with MyType \n";
      dc_ = dcConv_;
      dc_ = &(*dc_).operator += (dc);
    }
    else 
    {
      dc_ = const_cast<MyType *> (&dc);
      dcConv_ = const_cast<MyType *> (&dc);
    }
    return (*this);
  }
  
  //! Assignement operator
  template <class OpType>
  MyType & operator = (const OpType & dc)
  {
    //std::cout << "Store operator \n";
    dc_     = const_cast<OpType *> (&dc);
    dcConv_ = const_cast<OpType &> (dc).convert();
    return (*this);
  }

  //! Assignement operator
  MyType & operator = (const MyType & dc)
  {
    //std::cout << "No need to do this, use += \n";
    dc_     = const_cast<MyType *> (dc.dc_);
    dcConv_ = const_cast<MyType *> (dc.dcConv_);
    return (*this);
  }

  //! clear object list 
  virtual void clear() 
  {
    dc_ = 0;
    dcConv_ = 0;
  }
private: 
  MyType *dc_;
  MyType *dcConv_;
};


//! empty data collector 
template <class GridType>
class DummyDataCollector 
{
  typedef DummyDataCollector<GridType> MyType;
public:

  typedef std::pair < int * , int * > DataCollectorParamType;
  typedef LocalInterface<DataCollectorParamType> LocalInterfaceType;
  //! all adaptation operators have this method which adapts the
  //! corresponding grid and organizes the restriction prolongation process
  //! of the underlying function spaces
  void apply (int , int ) const 
  {
    std::cerr << "WARNING: apply: did nothing! \n";
  };
  
  //! Assignement operator
  template <class OpType>
  MyType & operator += (const OpType & dc)
  {
    return (*this);
  }
  
  //! Assignement operator
  template <class OpType>
  MyType & operator = (const OpType & dc)
  {
    return (*this);
  }
};

/** \brief The DataCollector is an example for a grid walk done while load
 * balancing moves entities from one processor to another. 
 * The Communicator or grid should call the inlineData (write Data to
 * ObjectStream) and the xtractData (read Data from Stream) and provide the
 * macro level Entity<codim =0> and the ObjectStream. This Operator then
 * does the hierarhic walk and calls its local pack operators which know
 * the discrete functions to pack to the stream. 
 */
template <class GridType, class LocalDataCollectImp>
class DataCollector
: public DataCollectorInterface<GridType, 
   typename GridObjectStreamOrDefault<GridType, DummyObjectStream > :: ObjectStreamType >  
, public ObjPointerStorage 
{  
  typedef typename GridObjectStreamOrDefault<GridType, DummyObjectStream > :: ObjectStreamType ObjectStreamType; 
  typedef typename GridType::template Codim<0>::Entity EntityType;
  typedef DataCollector<EntityType,LocalDataCollectImp> MyType;
  typedef DofManager<GridType> DofManagerType;

  typedef typename std::pair < ObjectStreamType * , const EntityType * > ParamType;
  typedef LocalInterface<ParamType> LocalInterfaceType;
  
  friend class DataCollectorInterface<GridType,ObjectStreamType>; 
  typedef DataCollectorInterface<GridType,ObjectStreamType> DataCollectorInterfaceType;
  
public:
  //! create DiscreteOperator with a LocalOperator 
  DataCollector (GridType & grid, DofManagerType & dm, LocalDataCollectImp & ldc, 
                 bool read , bool leaf , int numChildren = 8) 
    : grid_(grid) , dm_ ( dm ), ldc_ (ldc) 
    , rwType_((read) ? (readData) : writeData )
    , leaf_(leaf) 
    , numChildren_(numChildren) 
  {}

  //! Desctructor 
  virtual ~DataCollector () {}

  //! operator + (combine this operator) and return new Object 
  template <class LocalDataCollectType> 
  DataCollector<GridType,
  CombinedLocalDataCollect <LocalDataCollectImp,LocalDataCollectType> > & 
  operator + (const DataCollector<GridType,LocalDataCollectType> &op)
  {
    typedef DataCollector<GridType,LocalDataCollectType> CopyType;
    typedef CombinedLocalDataCollect <LocalDataCollectImp,LocalDataCollectType> COType;
     
    COType *newLDCOp = new COType ( ldc_  , const_cast<CopyType &> (op).getLocalOp() );
    typedef DataCollector <GridType, COType> OPType;
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp , (rwType_ == readData) );    

    // memorize this new generated object because is represents this
    // operator and is deleted if this operator is deleted
    saveObjPointer( dcOp , newLDCOp );
   
    return *dcOp;
  }

  //! oeprator += combine and return this Object 
  template <class LocalDataCollectType> 
  DataCollector<GridType,LocalInterface<ParamType> > & 
  operator += (const DataCollector<GridType,LocalDataCollectType> &op)
  {
    typedef DataCollector<GridType,LocalDataCollectType> CopyType;
    typedef LocalInterface<ParamType> COType;
     
    COType *newLDCOp = new COType ( ldc_ + op.getLocalOp() );
    typedef DataCollector <GridType, COType> OPType;
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp , (rwType_ == readData) );    

    // memorize this new generated object because is represents this
    // operator and is deleted if this operator is deleted
    saveObjPointer( dcOp , newLDCOp );
   
    return *dcOp;
  }
 
  //! operator += combine and return InterfaceType
  DataCollectorInterfaceType & 
  operator += (const DataCollectorInterfaceType &op)
  {
    //std::cout << "operator += with Interface Type \n";
    typedef LocalInterface<ParamType> COType;
    typedef DataCollector<GridType,COType> CopyType;
     
    COType *newLDCOp = new COType ( ldc_ + op.getLocalInterfaceOp() );
    typedef DataCollector<GridType, COType> OPType;
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp , (rwType_ == readData) , leaf_ );    

    // memorize this new generated object because is represents this
    // operator and is deleted if this operator is deleted
    saveObjPointer( dcOp , newLDCOp );
   
    return *dcOp;
  }
 
  //! return reference to loacl Operator 
  const LocalDataCollectImp & getLocalOp () const 
  {
    return ldc_;
  }
 
  //! return reference to loacl Operator 
  LocalDataCollectImp & getLocalOp ()  
  {
    return ldc_;
  }

  const LocalInterfaceType & getLocalInterfaceOp () const 
  {
    //std::cout << "getLocalInter \n";
    return ldc_;
  }
  
  LocalInterfaceType & getLocalInterfaceOp ()  
  {
    //std::cout << "getLocalInter \n";
    return ldc_;
  }

  //! apply, if this operator is in write status the inlineData is called
  //! else xtractData is called 
  void apply (ObjectStreamType & str, EntityType & en) const 
  {
    if(rwType_ == writeData) 
      inlineData(str,en);
    else 
      xtractData(str,en);
  }

  //! write all data of all entities blowe this Entity to the stream 
  void inlineData (ObjectStreamType & str, EntityType & en) const 
  {
    const int mxlvl = grid_.maxLevel();

    // read/write macro element
    inlineLocal(str,en);
    
    {
      typedef typename EntityType::HierarchicIterator HierarchicIteratorType;
      HierarchicIteratorType endit  = en.hend(mxlvl);
      for(HierarchicIteratorType it = en.hbegin(mxlvl); 
          it != endit; ++it )
      {
        inlineLocal(str, *it); 
      }
    }
  }

  //! read all data of all entities blowe this Entity from the stream 
  void xtractData (ObjectStreamType & str, EntityType & en) const 
  {
    const int mxlvl = grid_.maxLevel();

    // read/write macro element
    xtractLocal(str,en);
    
    {
      typedef typename EntityType::HierarchicIterator HierarchicIteratorType;
      HierarchicIteratorType endit  = en.hend(mxlvl);
      for(HierarchicIteratorType it = en.hbegin(mxlvl); 
          it != endit; ++it )
      {
        xtractLocal(str, *it); 
      }
    }
  }
  
private:
  enum ReadWriteType { readData , writeData };
  
  DataCollector<GridType,LocalInterface<ParamType> > * convert ()  
  {
    typedef LocalInterface<ParamType> COType;
     
    COType *newLDCOp = new COType ( ldc_ );
    typedef DataCollector <GridType, COType> OPType;
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp , (rwType_ == readData), leaf_ );    

    // memorize this new generated object because is represents this
    // operator and is deleted if this operator is deleted
    saveObjPointer( dcOp , newLDCOp );
   
    return dcOp;
  }
 
  // write data of entity 
  void inlineLocal(ObjectStreamType & str, EntityType& en) const 
  {
    assert( rwType_ == writeData );

    // if only leaf data is inlined then 
    // check whether en is leaf element 
    if( leaf_ && !en.isLeaf() ) return;

    ParamType p( &str , &en );
    // apply local operators 
    ldc_.apply( p );

    // remove entity from index sets 
    dm_.removeEntity( en );
  }
  
  // read data of entity 
  void xtractLocal(ObjectStreamType & str, EntityType& en) const 
  {
    assert( rwType_ == readData );
    
    // if only leaf data is inlined then 
    // check whether en is leaf element 
    if( leaf_ && !en.isLeaf() ) return;

    // insert entity into index sets 
    dm_.insertEntity( en );

    ParamType p( &str , &en );
    // apply local operators 
    ldc_.apply( p );
  }
  
  //! corresponding grid 
  mutable GridType & grid_;

  //! DofManager corresponding to grid
  mutable DofManagerType & dm_;
  
  //! Local Data Writer and Reader 
  mutable LocalDataCollectImp & ldc_;

  //! determines whether data is read or written
  const ReadWriteType rwType_;

  //! true if only leaf data is packed 
  const bool leaf_;
  
  // number of childs one element can have 
  const int numChildren_;
};


//***********************************************************************

    /** \brief ???
     * \todo Please doc me!
     */
template <class DiscreteFunctionType>
class DataInliner : 
public LocalInlinePlus < DataInliner< DiscreteFunctionType > , 
  typename std::pair < 
  typename GridObjectStreamOrDefault<typename DiscreteFunctionType::FunctionSpaceType::GridType, DummyObjectStream > :: ObjectStreamType * , 
      const typename DiscreteFunctionType::FunctionSpaceType::GridType::template Codim<0>::Entity * > > 
{
  typedef typename GridObjectStreamOrDefault<typename DiscreteFunctionType::FunctionSpaceType::GridType, DummyObjectStream > :: ObjectStreamType ObjectStreamType; 
  typedef LocalInlinePlus < DataInliner< DiscreteFunctionType > ,
    typename std::pair < ObjectStreamType * ,
          const typename  DiscreteFunctionType::FunctionSpaceType::GridType::template Codim<0>::Entity * > >  ChefType;
public:
  typedef typename DiscreteFunctionType::FunctionSpaceType::GridType::template Codim<0>::Entity EntityType;
  typedef typename ChefType::Traits::ParamType ParamType;

  typedef DataInliner<DiscreteFunctionType> MyType;

  typedef LocalInterface<ParamType> LocalInterfaceType;
    
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;

public:  
  DataInliner ( const DiscreteFunctionType & df ) 
    : df_ (df) 
  {}

  //! store data to stream  
  void apply ( ParamType & p ) const 
  {
    assert( p.first && p.second );
    this->apply( *p.first, *p.second );
  }

  //! store data to stream  
  void apply ( ObjectStreamType& str, const EntityType& en ) const 
  {
    assert( df_.space().indexSet().contains( en ) );
    
    const LocalFunctionType lf = df_.localFunction( en );
    const int numDofs = lf.numDofs();
    for(int l=0; l<numDofs; ++l)
    {
      str.write( lf[l] );
    }
  }

private:
  const DiscreteFunctionType & df_;
};


template <class DiscreteFunctionType>
class DataXtractor : 
public LocalInlinePlus < DataXtractor< DiscreteFunctionType > , 
  typename std::pair < 
  typename GridObjectStreamOrDefault<typename DiscreteFunctionType::FunctionSpaceType::GridType, DummyObjectStream > :: ObjectStreamType * , 
      const typename DiscreteFunctionType::FunctionSpaceType::GridType::template Codim<0>::Entity * > > 
{
  typedef typename GridObjectStreamOrDefault<typename DiscreteFunctionType::FunctionSpaceType::GridType, DummyObjectStream > :: ObjectStreamType ObjectStreamType;  
  typedef LocalInlinePlus < DataInliner< DiscreteFunctionType > ,
    typename std::pair < ObjectStreamType * ,
          const typename  DiscreteFunctionType::FunctionSpaceType::GridType::template Codim<0>::Entity * > >  ChefType;
public:
  typedef typename DiscreteFunctionType::FunctionSpaceType::GridType::template Codim<0>::Entity EntityType;
  typedef typename ChefType::Traits::ParamType ParamType;

  typedef DataXtractor<DiscreteFunctionType> MyType;

  typedef LocalInterface<ParamType> LocalInterfaceType;
    
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;
public:  
  DataXtractor ( DiscreteFunctionType & df ) 
    : df_ (df) 
    {}

  //! store data to stream  
  void apply ( ParamType & p ) const 
  {
    assert( p.first && p.second );
    this->apply( *p.first, *p.second );
  }

  //! store data to stream  
  void apply (ObjectStreamType & str, const EntityType & en ) const 
  {
    // make sure entity is contained in set 
    assert( df_.space().indexSet().contains( en ) );

    LocalFunctionType lf = df_.localFunction( en );
    const int numDofs = lf.numDofs();
    for(int l=0; l<numDofs; ++l)
    {
      str.read( lf[l] );
    }
  }

private:
  mutable DiscreteFunctionType & df_;
};

/** @} end documentation group */

}

#endif
