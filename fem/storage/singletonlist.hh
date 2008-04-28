#ifndef DUNE_SINGLETONLIST_HH
#define DUNE_SINGLETONLIST_HH

//- System includes 
#include <cassert>
#include <vector> 
#include <string>
#include <list>
#include <iostream>

namespace Dune
{

  template< class Key, class Object >
  struct DefaultSingletonFactory
  {
    static Object *createObject ( const Key &key )
    {
      return new Object( key );
    }

    static void deleteObject ( Object *object )
    {
      delete object;
    }
  };



  /** \class SingletonList
   *  \ingroup HelperClasses
   *  \brief Singleton list for key/object pairs
   *
   *  A singleton list guarantees that for any valid key at most one object is
   *  created.
   *
   *  \param  Key      type of keys
   *  \param  Object   type of objects
   *  \param  Factory  factory class creating objects from keys. The default
   *                   just passes the key to the object's constructor.
   */
  template< class Key, class Object,
            class Factory = DefaultSingletonFactory< Key, Object > >
  class SingletonList
  {
    typedef SingletonList< Key, Object, Factory > ThisType;

  public:
    typedef Key KeyType;
    typedef Object ObjectType;
    typedef Factory FactoryType;

    typedef std :: pair< ObjectType * , unsigned int * > ValueType;
    typedef std :: pair< KeyType, ValueType > ListObjType;

  private:
    typedef std :: list< ListObjType > ListType;
    typedef typename ListType :: iterator ListIteratorType;

    class SingletonListStorage;

  private:
    // prohibit creation
    SingletonList ();

    // prohibit copying
    SingletonList ( const ThisType & );

  public:
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
    inline static ObjectType &getObject( const KeyType &key )
    {
      ValueType objValue = getObjFromList( key );

      // if object exists, increase reference count and return it
      if( objValue.first )
      {
        ++( *(objValue.second) );
        return *(objValue.first);
      }

      // object does not exist. Create it with reference count of 1
      ObjectType *object = FactoryType :: createObject( key );
      assert( object );
      ValueType value( object, new unsigned int( 1 ) );
      ListObjType tmp( key, value ); 
      singletonList().push_back( tmp );
      return *object;
    } 

    //! decrease ref counter for this object, 
    //! if ref counter is zero, object is deleted 
    inline static void removeObject ( const ObjectType &object )
    {
      ListIteratorType end = singletonList().end();
      for( ListIteratorType it = singletonList().begin(); it != end; ++it )
      {
        ValueType value = (*it).second;
        if( value.first == &object )
        {
          eraseItem( it );
          return;
        }
      }
      std :: cerr << "Object could not be deleted, "
                  << "because it is not in the list anymore!" << std :: endl;
    }

    // return pair < Object * , refCounter *> 
    inline static ValueType getObjFromList( const KeyType &key )
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
    static void eraseItem( ListIteratorType &it )
    {
      ValueType value = (*it).second;
      unsigned int &refCount = *(value.second);

      assert( refCount > 0 );
      if( (--refCount) == 0 )
        deleteItem( it );
    }

  private:  
    static void deleteItem(ListIteratorType & it) 
    {
      ValueType val = (*it).second; 
      // remove from list
      singletonList().erase( it );
      // delete objects 
      FactoryType :: deleteObject( val.first );
      delete val.second;
    }
  }; // end SingletonList 



  
  template< class Key, class Object, class Factory >
  class SingletonList< Key, Object, Factory > :: SingletonListStorage
  {
    typedef SingletonListStorage ThisType;

  protected:
    ListType singletonList_; 

  public:  
    inline SingletonListStorage ()
    : singletonList_()
    {}
    
    inline ~SingletonListStorage ()
    {
      while( !singletonList().empty() )
        deleteItem( singletonList().begin() );
    }

    ListType &singletonList ()
    {
      return singletonList_;
    }

    void deleteItem ( const ListIteratorType &it )
    {
      ValueType val = (*it).second; 
      // remove from list
      singletonList().erase( it );
      // delete objects 
      FactoryType :: deleteObject( val.first );
      delete val.second;
    }
  };

} // end namespace Dune

#endif 
