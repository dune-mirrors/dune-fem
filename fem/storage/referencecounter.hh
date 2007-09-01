#ifndef DUNE_FEM_REFERENCECOUNTER_HH
#define DUNE_FEM_REFERENCECOUNTER_HH

#include <dune/common/bartonnackmanifcheck.hh>

namespace Dune
{

  template< class ReferenceCounterImp >
  class ReferenceCounterInterface
  {
  public:
    //! type of the implementation (Barton-Nackman)
    typedef ReferenceCounterImp ReferenceCounterType;

  private:
    typedef ReferenceCounterInterface< ReferenceCounterType > ThisType;

  public:
    typedef ThisType ReferenceCounterInterfaceType;

    template< class, class >
    friend class Conversion;

  public:
    inline ReferenceCounterInterface ()
    {
      typedef CompileTimeChecker
        < Conversion< ReferenceCounterType, ThisType > :: exists >
        __ReferenceCounter_Implementation_Must_Be_Derived_From_Interface__;
    }

    inline void addReference () const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().addReference() );
    }

    inline void deleteObject ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().deleteObject() );
    }

    inline void removeReference () const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().removeReference() );
    }

  protected:
    inline const ReferenceCounterType &asImp () const
    {
      return static_cast< const ReferenceCounterType & >( *this );
    }

    inline ReferenceCounterType &asImp ()
    {
      return static_cast< ReferenceCounterType & >( *this );
    }
  };



  template< class ReferenceCounterType >
  struct CheckReferenceCounterInterface
  {
    typedef ReferenceCounterInterface< ReferenceCounterType >
      ReferenceCounterInterfaceType;
    typedef CompileTimeChecker< Conversion< ReferenceCounterType,
                                            ReferenceCounterInterfaceType
                                          > :: exists >
      CheckerType;
  };



  template< class ReferenceCounterImp >
  class ReferenceCounterDefault
  : public ReferenceCounterInterface< ReferenceCounterImp >
  {
  public:
    typedef ReferenceCounterImp ReferenceCounterType;

  private:
    typedef ReferenceCounterDefault< ReferenceCounterType > ThisType;
    typedef ReferenceCounterInterface< ReferenceCounterType > BaseType;

    template< class, class >
    friend class Conversion;

  protected:
    using BaseType :: asImp;

  protected:
    mutable unsigned int refcount_;

  public:
    inline explicit ReferenceCounterDefault ( unsigned int refcount = 1 )
    : refcount_( refcount )
    {
    }

  private:
    // prohibit copying
    ReferenceCounterDefault ( const ThisType & );
    ThisType &operator= ( const ThisType &other );

  public:
    inline void addReference () const
    {
      ++refcount_;
    }

    inline void deleteObject ()
    {
      delete this;
    }

    inline void removeReference () const
    {
      assert( refcount_ > 0 );
      --refcount_;
      if( refcount_ == 0 )
        const_cast< ReferenceCounterType& >( asImp() ).deleteObject();
    }
  };



  template< class ObjectImp >
  class ObjectPointer
  {
  public:
    typedef ObjectImp ObjectType;

  private:
    typedef ObjectPointer< ObjectType > ThisType;

    typedef CheckReferenceCounterInterface< ObjectType > CheckObjectType;

  protected:
    ObjectType *object_;

  public:
    inline explicit ObjectPointer ( ObjectType *const object )
    : object_( object )
    {
      if( object_ != 0 )
        object_->addReference();
    }

    inline ObjectPointer ( const ThisType &other )
    : object_( other.object_ )
    {
      if( object_ != 0 )
        object_->addReference();
    }

    inline ~ObjectPointer ()
    {
      if( object_ != 0 )
        object_->removeReference();
    }

    inline ThisType &operator= ( const ThisType &other )
    {
      ObjectType *const newObject = other.object_;
      newObject->addReference();

      if( object_ != 0 )
        object_->removeReference();
      object_ = newObject;
    }

    ObjectType &operator* () const
    {
      return *object_;
    }
  };


 
  template< class ObjectImp >
  class ObjectReference
  {
  public:
    typedef ObjectImp ObjectType;

  private:
    typedef ObjectReference< ObjectType > ThisType;

    typedef CheckReferenceCounterInterface< ObjectType > CheckObjectType;

  protected:
    ObjectType &object_;

  public:
    inline explicit ObjectReference ( ObjectType &object )
    : object_( object )
    {
      object_.addReference();
    }

    inline ObjectReference ( const ThisType &other )
    : object_( other.object_ )
    {
      object_.addReference();
    }

    inline ~ObjectReference ()
    {
      object_.removeReference();
    }

  private:
    // prohibit copying
    ThisType &operator= ( const ThisType & );

  public:
    const ObjectType &operator* () const
    {
      return object_;
    }

    ObjectType &operator* ()
    {
      return object_;
    }
  };

}

#endif
