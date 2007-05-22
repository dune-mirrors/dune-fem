#ifndef DUNE_SINGLETONLIST_HH
#define DUNE_SINGLETONLIST_HH

//- System includes 
#include <cassert>
#include <vector> 
#include <string>
#include <list>

//- Dune includes 
#include <dune/common/dlist.hh>

namespace Dune { 

//! DefaultFactory that tells SingletonList how to create objects. 
//! Default is calling ObjectImp(*key) as constructor.
template <class KeyImp, class ObjectImp>
struct DefaultSingletonFactory  
{
  //! overload this method to create Objects with different constructors 
  static ObjectImp * createObject ( KeyImp key )
  {
    return new ObjectImp ( *key );
  }

  //! method for deleting pointer 
  static void deleteObject( ObjectImp * obj )
  {
    delete obj;
  }
};

//! Singleton list for key/object pairs that guarantees that for a given key only on object is created.
//! FactoryImp tells SingleList how to create objects and how to compare
//! keys. 
template <class KeyImp, class ObjectImp, class FactoryImp =
  DefaultSingletonFactory<KeyImp,ObjectImp> >
class SingletonList 
{
  typedef SingletonList<KeyImp,ObjectImp> ThisType;
    
  typedef KeyImp KeyType;
  typedef ObjectImp ObjectType;

  typedef std::pair < ObjectType * , size_t * > ValueType; 
  typedef std::pair < KeyType , ValueType > ListObjType;
  
  typedef DoubleLinkedList < ListObjType > ListType;
  typedef typename ListType::Iterator ListIteratorType;

  // destructor 
  class SingletonListStorage
  {
    ListType singletonList_; 
  public:  
    SingletonListStorage () {}
    
    ~SingletonListStorage () 
    { 
      ListIteratorType endit = singletonList().rend();
      ListIteratorType it = singletonList().rbegin(); 
      while( it != endit )
      {
        deleteItem( it );
        it = singletonList().rbegin();
      }
    }

    ListType& singletonList() { return singletonList_; }

    void deleteItem(ListIteratorType & it) 
    {
      ValueType val = (*it).second; 
      // remove from list
      singletonList().erase( it );
      // delete objects 
      FactoryImp::deleteObject( val.first );
      delete val.second;
    }
  };

  //! list that store pairs of key/object pointers 
  //! singleton list 
  inline static ListType & singletonList() 
  {
    static SingletonListStorage s; 
    //! list that store pairs of key/object pointers 
    return s.singletonList(); 
  }
  
public:  
  //! return reference to the object for given key. 
  //! If the object does not exist, then it is created first, otherwise the
  //! reference counter is increased. 
  //inline static ObjectType & getObject(KeyType key) 
  inline static ObjectType & getObject(const KeyType & key) 
  {
    // search list for dof manager 
    ValueType objVal = getObjFromList(key);

    // if not exists, create it, ref count is 1  
    if(!objVal.first)
    {
      ObjectType * obj = FactoryImp::createObject( key );
      assert( obj );
      // store pointer and ref count
      ValueType val ( obj , new size_t(1) );
      // store key and value 
      ListObjType tmp( key , val ); 
      singletonList().insert_after( singletonList().rbegin() , tmp ); 
      return *obj;
    }
    // if object exists, increase ref count 
    ++(*(objVal.second));
    return *(objVal.first);
  } 

  //! decrease ref counter for this object, 
  //! if ref counter is zero, object is deleted 
  inline static void removeObject (const ObjectType & obj) 
  {
    ListIteratorType endit = singletonList().end();
    for(ListIteratorType it = singletonList().begin(); it!=endit; ++it)
    {
      ValueType val = (*it).second; 
      if( val.first == & obj )
      {
        eraseItem(it);
        return;
      }
    }
    std::cerr << "Object could not be deleted, because it is not in the list anymore! \n";
  }

  // return pair < Object * , refCounter *> 
  inline static ValueType getObjFromList(const KeyType & key)
  {
    ListIteratorType endit = singletonList().end();
    for(ListIteratorType it = singletonList().begin(); it!=endit; ++it)
    {
      if( (*it).first == key )
      {
        return (*it).second; 
      }
    }
    return ValueType(0,0);
  }

protected:
  static void eraseItem(ListIteratorType & it) 
  {
    ValueType val = (*it).second; 
    assert( *(val.second) > 0 );
    // if only one reference left, remove object 
    if( *(val.second) == 1 )
    {
      // remove from list
      deleteItem( it );
    }
    else 
    {
      // decrease reference count 
      --(*(val.second));
    }
    return;
  }

private:  
  static void deleteItem(ListIteratorType & it) 
  {
    ValueType val = (*it).second; 
    // remove from list
    singletonList().erase( it );
    // delete objects 
    FactoryImp::deleteObject( val.first );
    delete val.second;
  }
  // constructor, ony this class is allowed to construct itself, prevent
  // from derivation 
  SingletonList () {};  
  SingletonList (const SingletonList &) {};  

}; // end SingletonList 

} // end namespace Dune 
#endif 
