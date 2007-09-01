#ifndef DUNE_OBJECTSTACK_HH
#define DUNE_OBJECTSTACK_HH

#include <dune/fem/misc/debug.hh>
#include <dune/fem/misc/referencecounter.hh>

namespace Dune
{

  template< class ObjectFactoryImp >
  class ObjectStack;



  template< class ObjectFactoryImp >
  class ObjectStackStorage
  : public ReferenceCounterDefault< ObjectStackStorage< ObjectFactoryImp > >
  {
  public:
    typedef ObjectFactoryImp ObjectFactoryType;

  private:
    typedef ObjectStackStorage< ObjectFactoryType > ThisType;
    typedef ReferenceCounterDefault< ThisType > BaseType;

    template< class, class >
    friend class Conversion;

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

    // next object on the stack
    ThisType *next_;

  protected:
    inline explicit ObjectStackStorage ( ObjectType *const obj,
                                         ObjectStackType &stack )
    : BaseType( 0 ),
      object_( obj ),
      stack_( stack )
    {
    }

  private:
    // prohibit copying
    ObjectStackStorage ( const ThisType & );

  public:
    inline ~ObjectStackStorage ()
    {
      delete object_;
    }
    
  private:
    // prohivit assignment
    ThisType &operator= ( const ThisType & );

  public:
    inline operator const ObjectType& () const
    {
      return *object_;
    }

    inline operator ObjectType& ()
    {
      return *object_;
    }

    inline void deleteObject ()
    {
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
    typedef ObjectReference< StorageType > ObjectReferenceType;
    
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
    ObjectStack ( const ThisType & );

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
    ThisType &operator= ( const ThisType & );

  public:
    //! get an object pointer to a storage object
    inline ObjectReferenceType getObject ()
    {
      return ObjectReferenceType( pop() );
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
        ptr = new StorageType( factory().newObject(), *this );
      }
      return *ptr;
    }

  protected:
    inline const ObjectFactoryType &factory() const
    {
      return factory_;
    }
  };

} // end namespace Dune

#endif
