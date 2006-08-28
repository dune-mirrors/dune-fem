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

//! DofManagerFactory guarantees that only one instance of a dofmanager 
//! per grid is generated. If getDofManager is called with a grid for which
//! already a mamager exists, then the reference to this manager is returned. 
template <class KeyImp, class ObjectImp>
class SingletonList 
{
  typedef SingletonList<KeyImp,ObjectImp> ThisType;
    
  typedef KeyImp KeyType;
  typedef ObjectImp ObjectType;

  typedef std::pair < ObjectType * , size_t * > ValueType; 
  typedef std::pair < const KeyType * , ValueType > ListObjType;
  
  typedef DoubleLinkedList < ListObjType > ListType;
  typedef typename ListType::Iterator ListIteratorType;

  //! list that store pairs of key/object pointers 
  //! singleton list 
  inline static ListType & singletonList() 
  {
    //! list that store pairs of key/object pointers 
    static ListType singletonList_; 
    return singletonList_; 
  }
  
public:  
  //! return reference to the DofManager for the given grid. 
  //! If the object does not exist, then it is created first, otherwise the
  //! reference counter is increased 
  inline static ObjectType & getObject(const KeyType & key) 
  {
    // search list for dof manager 
    ValueType objVal = getObjFromList(key);

    // if not exists, create it, ref count is 1  
    if(!objVal.first)
    {
      ObjectType * obj = new ObjectType ( key );
      assert( obj );
      // store pointer and ref count
      ValueType val ( obj , new size_t(1) );
      // store key and value 
      ListObjType tmp( &key , val ); 
      singletonList().insert_after( singletonList().rbegin() , tmp ); 
      return *obj;
    }
    // if object exists, increase ref count 
    ++(*(objVal.second));
    return *(objVal.first);
  } 

  //! delete the dof manager that belong to the given grid 
  inline static void removeObject (const ObjectType & obj) 
  {
    ListIteratorType endit = singletonList().end();
    for(ListIteratorType it = singletonList().begin(); it!=endit; ++it)
    {
      ValueType val = (*it).second; 
      if( val.first == (& obj))
      {
        eraseItem(it);
        return;
      }
    }
    std::cerr << "Object could not be deleted, because it is not in the list anymore! \n";
  }

  // return pointer to dof manager for given grid 
  inline static ValueType getObjFromList(const KeyType & key)
  {
    ListIteratorType endit = singletonList().end();
    for(ListIteratorType it = singletonList().begin(); it!=endit; ++it)
    {
      if( (*it).first == & key )
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
    //std::cout << (*(val.second)) << " ref count is \n";
    if( *(val.second) == 1 )
    {
      singletonList().erase( it );
      //std::cout << "Deleting obj = " << val.first << "\n";
      delete val.first;
      delete val.second;
    }
    else 
    {
      //std::cout << "descreasing refcount for obj = " << val.first << "\n";
      --(*(val.second));
      //std::cout << (*(val.second)) << " ref count is \n";
    }
    return;
  }

private:  
  // constructor, ony this class is allowed to construct itself, prevent
  // from derivation 
  SingletonList () {};  
  SingletonList (const SingletonList &) {};  

  // destructor 
  ~SingletonList () 
  { 
  }

}; // end SingletonList 

} // end namespace Dune 
#endif 
