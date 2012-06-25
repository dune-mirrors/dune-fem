#ifndef DUNE_DATACOLLECTOR_HH
#define DUNE_DATACOLLECTOR_HH

//-System includes
#include <cassert>
#include <cstdlib>

#include <vector>
#include <utility>
#include <iostream>

//-Dune includes#
#include <dune/common/version.hh>

//- local includes 
#include <dune/fem/operator/common/objpointer.hh>

namespace Dune{

namespace Fem {

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


// define read or write stream 
struct DataCollectorTraits 
{
  enum ReadWriteType { readData , writeData };
};

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

  template <class PT> 
  struct ObjectStreamExtractor 
  {
    typedef PT ObjectStreamType ;
  };

  template < class T1 , class T2 > 
  struct ObjectStreamExtractor< std::pair< T1* , const T2* > > 
  {
    typedef T1 ObjectStreamType ;
  };

  typedef typename ObjectStreamExtractor< ParamT > :: ObjectStreamType ObjectStreamType;

protected:
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
  static void copyList (ListType & vec, const MyType & op )
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
  //ReadWriteType rwType_ ; 
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
    this->saveObjPointer( combo );
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
public:  
  typedef ObjectStreamImp  ObjectStreamType;


  typedef typename GridType::template Codim<0>::Entity EntityType;

  typedef DataCollectorInterface<GridType, ObjectStreamImp> MyType;
  typedef std::pair < ObjectStreamType * , const EntityType * >  DataCollectorParamType;
  
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
  virtual void apply (ObjectStreamType &str, const EntityType & entity ) const 
  {
    //std::cout << "apply on interface class \n";
    if(dc_) (*dc_).apply(str, entity );  
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
template <class GridType, 
          class LocalDataCollectImp > 
class DataCollector
: public DataCollectorInterface< GridType, typename LocalDataCollectImp :: ObjectStreamType >  
, public ObjPointerStorage 
{  
public:  
  typedef typename LocalDataCollectImp :: ObjectStreamType ObjectStreamType; 
  typedef typename GridType::template Codim<0>::Entity EntityType;

protected:  
  typedef DataCollectorInterface< GridType, ObjectStreamType > BaseType;
  typedef typename DataCollectorTraits :: ReadWriteType ReadWriteType;

  typedef DataCollector<EntityType,LocalDataCollectImp> MyType;
  typedef DofManager<GridType> DofManagerType;

  typedef typename std::pair < ObjectStreamType * , const EntityType * > ParamType;
  typedef LocalInterface<ParamType> LocalInterfaceType;
  
  friend class DataCollectorInterface<GridType, ObjectStreamType>; 
  typedef DataCollectorInterface<GridType, ObjectStreamType> DataCollectorInterfaceType;
  
public:
  //! create DiscreteOperator with a LocalOperator 
  DataCollector (GridType & grid, 
                 DofManagerType & dm, 
                 LocalDataCollectImp & ldc, 
                 const ReadWriteType rwType,
                 int numChildren = 8) 
    : grid_(grid) , dm_ ( dm ), ldc_ (ldc) 
    , rwType_( rwType )
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
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp, rwType_ );    

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
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp, rwType_ );    

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
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp, rwType_ );  

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

  //! return true if data collector is writing data instead of reading 
  bool writeData() const  { return rwType_ == DataCollectorTraits :: writeData ; }

  //! apply, if this operator is in write status the inlineData is called
  //! else xtractData is called 
  void apply ( ObjectStreamType & str, const EntityType & entity ) const 
  {
    if( writeData() ) 
      inlineData(str, entity );
    else 
      xtractData(str, entity );
  }

  //! write all data of all entities blowe this Entity to the stream 
  void inlineData (ObjectStreamType & str, const EntityType & entity ) const 
  {
    const int mxlvl = grid_.maxLevel();

    // read/write macro element
    inlineLocal(str, entity );
    
    {
      typedef typename EntityType::HierarchicIterator HierarchicIteratorType;
      const HierarchicIteratorType endit  = entity.hend( mxlvl );
      for(HierarchicIteratorType it = entity.hbegin( mxlvl ); 
          it != endit; ++it )
      {
        inlineLocal(str, *it); 
      }
    }
  }

  //! read all data of all entities blowe this Entity from the stream 
  void xtractData (ObjectStreamType & str, const EntityType & entity ) const 
  {
    const int mxlvl = grid_.maxLevel();

    // read/write macro element
    xtractLocal(str, entity );
    
    {
      typedef typename EntityType::HierarchicIterator HierarchicIteratorType;
      const HierarchicIteratorType endit  = entity.hend( mxlvl );
      for(HierarchicIteratorType it = entity.hbegin( mxlvl ); 
          it != endit; ++it )
      {
        xtractLocal(str, *it); 
      }
    }
  }
  
private:
  DataCollector<GridType,LocalInterface<ParamType> > * convert ()  
  {
    typedef LocalInterface<ParamType> COType;
     
    COType *newLDCOp = new COType ( ldc_ );
    typedef DataCollector <GridType, COType> OPType;
   
    OPType *dcOp = new OPType ( grid_ , dm_ , *newLDCOp, rwType_ );    

    // memorize this new generated object because is represents this
    // operator and is deleted if this operator is deleted
    saveObjPointer( dcOp , newLDCOp );
   
    return dcOp;
  }
 
  // write data of entity 
  void inlineLocal(ObjectStreamType & str, const EntityType& entity ) const 
  {
    assert( writeData() );

    ParamType p( &str , &entity );
    // apply local operators 
    ldc_.apply( p );
  }
  
  // read data of entity 
  void xtractLocal(ObjectStreamType & str, const EntityType& entity ) const 
  {
    assert( ! writeData() );
    
    ParamType p( &str , &entity );
    // apply local operators 
    ldc_.apply( p );
  }
  
  //! corresponding grid 
  GridType &grid_;

  //! DofManager corresponding to grid
  DofManagerType &dm_;
  
  //! Local Data Writer and Reader 
  LocalDataCollectImp &ldc_;

  //! determines whether data is read or written
  const ReadWriteType rwType_;

  // number of childs one element can have 
  const int numChildren_;
};


//***********************************************************************
template< class DiscreteFunctionType >
struct LocalDataInlinerTraits
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridType        GridType;
  typedef typename DiscreteFunctionSpaceType::EntityType      EntityType;
  typedef typename GridType :: template Codim< 0 > :: Entity  GridEntityType;

  typedef DofManager < GridType > DofManagerType ;
  typedef typename DofManagerType :: InlineStreamType ObjectStreamType;

  typedef std::pair< ObjectStreamType *, const GridEntityType * > ParamType;
  typedef LocalInterface< ParamType > LocalInterfaceType;
};


/** \brief ???
 * \todo Please doc me!
 */
template< class DiscreteFunctionType, 
          class ContainsCheck >
class LocalDataInliner
: public LocalInlinePlus< LocalDataInliner< DiscreteFunctionType, ContainsCheck >,
                          typename LocalDataInlinerTraits< DiscreteFunctionType >::ParamType >
{
  typedef LocalDataInliner< DiscreteFunctionType, ContainsCheck > ThisType;
  typedef LocalInlinePlus< ThisType, typename LocalDataInlinerTraits< DiscreteFunctionType >::ParamType > BaseType;

public:
  typedef LocalDataInlinerTraits< DiscreteFunctionType > Traits;
  typedef typename Traits::ObjectStreamType ObjectStreamType;

  typedef typename Traits :: DofManagerType  DofManagerType;

  typedef typename Traits::EntityType     EntityType;
  typedef typename Traits::GridEntityType GridEntityType;
  typedef typename Traits::ParamType ParamType;

  typedef LocalInterface<ParamType> LocalInterfaceType;
    
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;

  //! constructor 
  LocalDataInliner ( const DiscreteFunctionType & df, 
                     const ContainsCheck& containsCheck ) 
    : df_ (df),
      dm_( DofManagerType :: instance( df.space().grid() ) ),
      containsCheck_( containsCheck )
  {}

  //! copy constructor 
  LocalDataInliner ( const LocalDataInliner & other  ) 
    : df_ (other.df_),
      dm_( other.dm_ ),
      containsCheck_( other.containsCheck_ )
  {}

  //! store data to stream  
  void apply ( ParamType & p ) const 
  {
    assert( p.first && p.second );
    const EntityType& entity = df_.space().gridPart().convert( *p.second );
    inlineData( *p.first, entity, *p.second );
  }

  typename DataCollectorTraits :: ReadWriteType
  readWriteInfo() const { return DataCollectorTraits :: writeData ; }
protected:
  //! store data to stream  
  void inlineData ( ObjectStreamType& str, 
                    const EntityType& entity, 
                    const GridEntityType& gridEntity ) const 
  {
    if( ! containsCheck_.contains ( entity ) ) return ;

    assert( df_.space().indexSet().contains( entity ) );
    
    const LocalFunctionType lf = df_.localFunction( entity );
    const int numDofs = lf.numDofs();
    for(int l=0; l<numDofs; ++l)
    {
      str.write( lf[l] );
    }

    // remove entity from index sets 
    dm_.removeEntity( gridEntity );
  }

protected:
  const DiscreteFunctionType & df_;
  DofManagerType& dm_ ;
  const ContainsCheck containsCheck_;
};


template< class DiscreteFunctionType >
struct LocalDataXtractorTraits
{
  typedef typename DiscreteFunctionType::DiscreteFunctionSpaceType
    DiscreteFunctionSpaceType;
  typedef typename DiscreteFunctionSpaceType::GridType        GridType;
  typedef typename DiscreteFunctionSpaceType::EntityType      EntityType;
  typedef typename GridType :: template Codim< 0 > :: Entity  GridEntityType;

  typedef DofManager < GridType > DofManagerType ;
  typedef typename DofManagerType :: XtractStreamType ObjectStreamType;

  typedef std::pair< ObjectStreamType *, const GridEntityType * > ParamType;
  typedef LocalInterface< ParamType > LocalInterfaceType;
};


template< class DiscreteFunctionType,
          class ContainsCheck >
class LocalDataXtractor
: public LocalInlinePlus< LocalDataXtractor< DiscreteFunctionType, ContainsCheck >,
                          typename LocalDataXtractorTraits< DiscreteFunctionType >::ParamType >
{
  typedef LocalDataXtractor< DiscreteFunctionType, ContainsCheck > ThisType;
  typedef LocalInlinePlus< ThisType, typename LocalDataXtractorTraits< DiscreteFunctionType >::ParamType > BaseType;

public:
  typedef LocalDataXtractorTraits< DiscreteFunctionType > Traits;
  typedef typename Traits::ObjectStreamType ObjectStreamType;

  typedef typename Traits :: DofManagerType  DofManagerType;

  typedef typename Traits::EntityType     EntityType;
  typedef typename Traits::GridEntityType GridEntityType;
  typedef typename Traits::ParamType ParamType;

  typedef LocalInterface<ParamType> LocalInterfaceType;
    
  typedef typename DiscreteFunctionType::LocalFunctionType LocalFunctionType;

  typedef typename DiscreteFunctionType::RangeFieldType RangeFieldType;
  typedef typename DiscreteFunctionType::DomainType DomainType;

  //! constructor 
  LocalDataXtractor ( DiscreteFunctionType & df, 
                      const ContainsCheck& containsCheck ) 
    : df_ (df),
      dm_( DofManagerType :: instance( df.space().grid() ) ),
      containsCheck_( containsCheck )
    {}

  //! copy constructor 
  LocalDataXtractor ( const LocalDataXtractor & other ) 
    : df_( other.df_ ),
      dm_( other.dm_ ),
      containsCheck_( other.containsCheck_ )
  {}

  //! store data to stream  
  void apply ( ParamType & p ) const 
  {
    assert( p.first && p.second );
    const EntityType& entity = df_.space().gridPart().convert( *p.second );
    xtractData( *p.first, entity, *p.second );
  }
  
  typename DataCollectorTraits :: ReadWriteType
  readWriteInfo() const { return DataCollectorTraits :: readData ; }
protected:
  //! store data to stream  
  void xtractData (ObjectStreamType & str, 
                   const EntityType& entity,
                   const GridEntityType& gridEntity ) const 
  {
    if( ! containsCheck_.contains ( entity ) ) return ;

    // insert entity into index sets 
    dm_.insertEntity( gridEntity );

    // make sure entity is contained in set 
    assert( df_.space().indexSet().contains( entity ) );

    LocalFunctionType lf = df_.localFunction( entity );
    const int numDofs = lf.numDofs();
    for(int l=0; l<numDofs; ++l)
    {
      str.read( lf[l] );
    }
  }

protected:
  DiscreteFunctionType &df_;
  DofManagerType &dm_;
  const ContainsCheck containsCheck_;
};

/** @} end documentation group */

} // namespace Fem 

} // namespace Dune 

#endif
