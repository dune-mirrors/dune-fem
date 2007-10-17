#ifndef DUNE_OBJECTSTACK_HH
#define DUNE_OBJECTSTACK_HH

#include <dune/fem/misc/debug.hh>
#include <dune/fem/storage/referencecounter.hh>

namespace Dune
{

  template< class ObjectFactoryImp >
  class ObjectStack;

  template< class ObjectFactoryImp >
  class ObjectStackEntry;


  
  template< class ObjectFactoryImp >
  struct ObjectStackEntryTraits
  {
    typedef ObjectStackEntry< ObjectFactoryImp > ReferenceCounterType;

    typedef typename ObjectFactoryImp :: ObjectType ObjectType;
  };

  

  template< class ObjectFactoryImp >
  class ObjectStackEntry
  : public ReferenceCounterDefault< ObjectStackEntryTraits< ObjectFactoryImp > >
  {
  public:
    typedef ObjectFactoryImp ObjectFactoryType;

    typedef ObjectStackEntryTraits< ObjectFactoryType > Traits;

  private:
    typedef ObjectStackEntry< ObjectFactoryType > ThisType;
    typedef ReferenceCounterDefault< Traits > BaseType;

    template< class, class >
    friend class Conversion;

    friend class ObjectStack< ObjectFactoryType >;

  protected:
    typedef ObjectStack< ObjectFactoryType > ObjectStackType;

  public:
    typedef typename ObjectFactoryType :: ObjectType ObjectType;

  protected:
    // reference to the stack
    ObjectStackType &stack_;

    // pointer to the actual object
    ObjectType *const object_;

    // next object on the stack
    ThisType *next_;

  protected:
    inline explicit ObjectStackEntry ( ObjectStackType &stack )
    : BaseType( 0 ),
      stack_( stack ),
      object_( stack_.factory().newObject() )
    {
    }

  private:
    // prohibit copying
    ObjectStackEntry ( const ThisType & );

  public:
    inline ~ObjectStackEntry ()
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
      stack_.push( this );
    }

    inline const ObjectType &getObject () const
    {
      return *object_;
    }
    
    inline ObjectType &getObject ()
    {
      return *object_;
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

    friend class ObjectStackEntry< ObjectFactoryType >;

  public:
    //! type of the stored objects
    typedef typename ObjectFactoryType :: ObjectType ObjectType;

    //! type of the storage objects
    typedef ObjectStackEntry< ObjectFactoryType > StackEntryType;

    //! type of object pointers
    typedef ObjectPointer< StackEntryType > PointerType;
    
  protected:
    const ObjectFactoryType &factory_;
    StackEntryType *top_;

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
        StackEntryType *obj = top_;
        top_ = top_->next_;
        delete obj;
      }
    }
    
  private:
    // Disallow copying
    ThisType &operator= ( const ThisType & );

  public:
    //! get an object pointer to a storage object
    inline PointerType getObject ()
    {
      return PointerType( pop() );
    }

  protected:
    inline const ObjectFactoryType &factory() const
    {
      return factory_;
    }

    // push storage object to the stack
    inline void push ( StackEntryType *obj )
    {
      --numIssuedObjects_;
      obj->next_ = top_;
      top_ = obj;
    }

    // pop a storage object from the stack
    inline StackEntryType *pop ()
    {
      ++numIssuedObjects_;

      StackEntryType *ptr = top_;
      if( ptr != 0 )
        top_ = top_->next_;
      else
        ptr = new StackEntryType( *this );
      return ptr;
    }
  };

} // end namespace Dune

#endif
