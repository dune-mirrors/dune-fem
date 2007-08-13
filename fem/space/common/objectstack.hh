#ifndef DUNE_OBJECTSTACK_HH
#define DUNE_OBJECTSTACK_HH

#include <vector>
#include <stack>

#include <dune/fem/misc/debug.hh>

namespace Dune
{

  template< class StorageImp >
  class ObjectPointer
  {
  public:
    typedef StorageImp StorageType;

  private:
    typedef ObjectPointer< StorageType > ThisType;

  public:
    typedef typename StorageType :: ObjectType ObjectType;

  protected:
    StorageType &storage_;

  public:
    inline explicit ObjectPointer ( StorageType &storage )
    : storage_( storage )
    {
      storage_.addReference();
    }

    inline ObjectPointer ( const ThisType &other )
    : storage_( other.storage_ )
    {
      storage_.addReference();
    }

    inline ~ObjectPointer ()
    {
      storage_.removeReference();
    }
    
  private:
    // Disallow assigning
    inline ThisType &operator= ( const ThisType &other );

  public:
    //! derefenence pointer
    ObjectType &operator* ()
    {
      return storage_.object();
    }

    //! dereference pointer
    const ObjectType &operator* () const
    {
      return storage_.object();
    }
  };



  template< class ObjectFactoryImp >
  class ObjectStack;



  template< class ObjectFactoryImp >
  class ObjectStackStorage
  {
  public:
    typedef ObjectFactoryImp ObjectFactoryType;

  private:
    typedef ObjectStackStorage< ObjectFactoryType > ThisType;

    friend class ObjectStack< ObjectFactoryType >;

  protected:
    typedef ObjectStack< ObjectFactoryType > ObjectStackType;

  public:
    typedef typename ObjectFactoryType :: ObjectType ObjectType;

  protected:
    // pointer to the actual object
    ObjectType *const object_;

    // reference to the stack
    ObjectStackType &stack_;

    // reference counter to this object
    unsigned int refcount_;
    
    // next object on the stack
    ThisType *next_;

  protected:
    inline explicit ObjectStackStorage ( ObjectType *const obj,
                                         ObjectStackType &stack )
    : object_( obj ),
      stack_( stack ),
      refcount_( 0 )
    {
    }

  private:
    // Disallow copying
    inline ObjectStackStorage ( const ThisType &other );

  public:
    inline ~ObjectStackStorage ()
    {
      assert( refCount_ == 0 );
      delete object_;
    }
    
  private:
    // Disallow copying
    inline ThisType &operator= ( const ThisType &other );

  public:
    inline ObjectType &object ()
    {
      return *object_;
    }

    inline void addReference ()
    {
      ++refcount_;
    }

    inline void removeReference ()
    {
      assert( refcount_ > 0 );
      --refcount_;
      if( refcount_ == 0 )
        stack_.push( *this );
    }
  };



  //! \ingroup HelperClasses
  //! Stores pointers to a given class in a stack
  //! used for local functions and for basefunctionsets
  template< class ObjectFactoryImp >
  class ObjectStack 
  {
    typedef ObjectFactoryImp ObjectFactoryType;

  private:
    typedef ObjectStack< ObjectFactoryType > ThisType;

  public:
    //! type of the stored objects
    typedef typename ObjectFactoryType :: ObjectType ObjectType;

    //! type of the storage objects
    typedef ObjectStackStorage< ObjectFactoryType > StorageType;

    //! type of object pointers
    typedef ObjectPointer< StorageType > ObjectPointerType;
    
  protected:
    const ObjectFactoryType &factory_;
    StorageType *top_;

    DebugCounter<> numIssuedObjects_;

  public:
    //! constructor 
    ObjectStack ( const ObjectFactoryType &factory )
    : factory_( factory ),
      top_( 0 )
    {
    }

  private:
    // Disallow copying
    ObjectStack ( const ThisType &other );

  public:
    //! delete all objects on stack 
    ~ObjectStack ()
    {
      assert( numIssuedObjects_ == 0 );

      while ( top_ != 0 )
      {
        StorageType *obj = top_;
        top_ = top_->next_;
        delete obj;
      }
    }
    
  private:
    // Disallow copying
    ThisType &operator= ( const ThisType &other );

  public:
    //! get an object pointer to a storage object
    inline ObjectPointerType getObject ()
    {
      return ObjectPointerType( pop() );
    }

    //! push storage object to the stack
    inline void push ( StorageType &obj )
    {
      --numIssuedObjects_;
      obj.next_ = top_;
      top_ = &obj;
    }

    //! pop a storage object from the stack
    inline StorageType &pop ()
    {
      ++numIssuedObjects_;

      StorageType *ptr = top_;
      if( ptr != 0 )
        top_ = top_->next_;
      else {
        ptr = new StorageType( factory_.newObject(), *this );
      }
      return *ptr;
    }
  };

} // end namespace Dune

#endif
