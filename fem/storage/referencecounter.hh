#ifndef DUNE_FEM_REFERENCECOUNTER_HH
#define DUNE_FEM_REFERENCECOUNTER_HH

#include <cassert>

#include <dune/fem/misc/bartonnackmaninterface.hh>

namespace Dune
{

  /** \class ReferenceCounterInterface
   *  \brief interface for objects capable of reference counting
   *
   *  In many situations it has to be decided, whether an object is still
   *  required or must be freed. A simple approach to this decision is to
   *  count all references to an object. When the last reference is released,
   *  the object is freed.
   *  
   *  The ReferenceCounterInterface provides methods, an object needs to
   *  support reference counting. Classes implementing this interface can be
   *  used with ObjectReference and ObjectPointer. Both will implicitly call
   *  addReference when they start pointing to an object and removeReference
   *  when they stop pointing to it.
   *
   *  \note Reference counting is not a good mehtod for garbage collection,
   *        since it does not detect cycles. So be careful when using reference
   *        counters!
   */
  template< class TraitsImp >
  class ReferenceCounterInterface
  : public BartonNackmanInterface< ReferenceCounterInterface< TraitsImp >,
                                   typename TraitsImp :: ReferenceCounterType >
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

    //! type of the implementation (Barton-Nackman)
    typedef typename Traits :: ReferenceCounterType ReferenceCounterType;

  private:
    typedef ReferenceCounterInterface< Traits > ThisType;
    typedef BartonNackmanInterface< ThisType, ReferenceCounterType > BaseType;

    template< class, class >
    friend class Conversion;

  public:
    //! type of the reference counter interface
    typedef ThisType ReferenceCounterInterfaceType;

    //! type of the object, this is a reference counter for
    typedef typename Traits :: ObjectType ObjectType;

  protected:
    using BaseType :: asImp;

  public:
    /** \brief add a reference to this object
     *
     *  This method should be called whenever a permanent reference to this
     *  object is established (a pointer counts as a reference, too)
     *
     *  \note This method is declared const, since we want to be able to add
     *        references to const objects, too. Hence, the reference coutner
     *        must be declared mutable.
     */
    inline void addReference () const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().addReference() );
    }

    /** \brief delete to object
     *
     *  This method is used by the default implementation to free the object
     *  when the reference counter becomes zero.
     *
     *  \note For the implementor's convenience, this method is not declared
     *        const. The caller (usually the default implementation's
     *        removeReference method) is responsible for performing a necessary
     *        const_cast.
     */
    inline void deleteObject ()
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().deleteObject() );
    }

    /** \brief access the real object (const version)
     *
     *  Sometimes the reference counter just wraps the object for which it
     *  counts references. In such cases, this method returns the wrapped
     *  object. Otherwise, the object itself may be returned.
     *
     *  \returns a constant reference to the real object
     */
    inline const ObjectType &getObject () const
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().getObject() );
      return asImp().getObject();
    }

    /** \brief access the real object (non-const version)
     *
     *  Sometimes the reference counter just wraps the object for which it
     *  counts references. In such cases, this method returns the wrapped
     *  object. Otherwise, the object itself may be returned.
     *
     *  \returns a reference to the real object
     */
    inline ObjectType &getObject ()
    {
      CHECK_INTERFACE_IMPLEMENTATION( asImp().getObject() );
      return asImp().getObject();
    }

    /** \brief remove a reference to this object
     *
     *  This method should be called whenever a previously added reference
     *  to this object is released.
     *
     *  \note This method is declared const, since we want to be able to remove
     *        references to const objects, too. Hence, the reference coutner
     *        must be declared mutable.
     */
    inline void removeReference () const
    {
      CHECK_AND_CALL_INTERFACE_IMPLEMENTATION( asImp().removeReference() );
    }
  };



  template< class ReferenceCounterType >
  struct CheckReferenceCounterInterface
  {
    typedef ReferenceCounterInterface< typename ReferenceCounterType :: Traits >
      ReferenceCounterInterfaceType;
    typedef CompileTimeChecker< Conversion< ReferenceCounterType,
                                            ReferenceCounterInterfaceType
                                          > :: exists >
      CheckerType;
  };



  /** \class ReferenceCounterDefault
   *  \brief default implementation of ReferenceCounterInterface
   *
   *  Reference counting usually uses a class variable to hold the number of
   *  references to the instance. To make reference counting simple to use,
   *  this mechanism is default implemented.
   *
   *  \note The deleteObject method is used to delete the object when the
   *        reference counter reaches zero. To perform any action other than
   *        deleting the object, simply overwrite this method.
   */
  template< class TraitsImp >
  class ReferenceCounterDefault
  : public ReferenceCounterInterface< TraitsImp >
  {
  public:
    //! type of the traits
    typedef TraitsImp Traits;

  private:
    typedef ReferenceCounterDefault< Traits > ThisType;
    typedef ReferenceCounterInterface< Traits > BaseType;

    template< class, class >
    friend class Conversion;

  public:
    //! type of the implementation (Barton-Nackman)
    typedef typename Traits :: ReferenceCounterType ReferenceCounterType;

  protected:
    using BaseType :: asImp;

  protected:
    mutable unsigned int refcount_;

  public:
    /** \brief constructor initializing the reference counter
     *
     *  \note Since we want to be able to count referenced to stack objects,
     *        we initialize the reference counter to 1 by default. This means
     *        that after calling new we already have one reference to the
     *        object. In such cases, just call removeReference directly when
     *        the reference is no longer used (instead of simply deleting it).
     *
     *  \param[in]  refcount  start value for the reference counter; this
     *              value defalts to 1
     */
    inline explicit ReferenceCounterDefault ( unsigned int refcount = 1 )
    : refcount_( refcount )
    {
    }

  private:
    // prohibit copying
    ReferenceCounterDefault ( const ThisType & );
    ThisType &operator= ( const ThisType &other );

  public:
    /** \copydoc Dune :: ReferenceCounterInterface :: addReference */
    inline void addReference () const
    {
      ++refcount_;
    }

    /** \copydoc Dune :: ReferenceCounterInterface :: deleteObject */
    inline void deleteObject ()
    {
      delete this;
    }

    /** \copydoc Dune :: ReferenceCounterInterface :: removeReference */
    inline void removeReference () const
    {
      assert( refcount_ > 0 );
      --refcount_;
      if( refcount_ == 0 )
        const_cast< ReferenceCounterType& >( asImp() ).deleteObject();
    }
  };



  /** \class ObjectPointer
   *  \brief models a pointer to a reference countable object
   *
   *  ObjectPointer tries to behave like a normal pointer to an object
   *  implementing the ReferenceCounterInterface. Internally, however, it
   *  calls the object's addReference and removeReference methods whenever the
   *  pointer is created, assign or deleted.
   */
  template< class ReferenceCounterImp >
  class ObjectPointer
  {
  public:
    //! type of the object, this pointer points to
    typedef ReferenceCounterImp ReferenceCounterType;

  private:
    typedef ObjectPointer< ReferenceCounterType > ThisType;

    typedef CheckReferenceCounterInterface< ReferenceCounterType >
      CheckReferenceCounterType;

  public:
    typedef typename ReferenceCounterType :: ObjectType ObjectType;

  protected:
    ReferenceCounterType *object_;

  public:
    /** \brief initialize a pointer (with a standard C++ pointer)
     *
     *  \param[in]  object  C++ pointer to initialize this pointer with; the
     *                      default value is 0
     */
    inline explicit ObjectPointer ( ReferenceCounterType *const object = 0 )
    : object_( object )
    {
      if( object_ != 0 )
        object_->addReference();
    }

    /** \brief copy constructor
     *
     *  Copying an ObjectPointer will also increase the reference counter of
     *  the object pointed to.
     *
     *  \param[in]  other  pointer to assign to this one
     */
    inline ObjectPointer ( const ThisType &other )
    : object_( other.object_ )
    {
      if( object_ != 0 )
        object_->addReference();
    }

    /** \brief destructor
     *
     *  When the pointer is deleted, the reference counter of the object pointed
     *  to is automatically decreased.
     */
    inline ~ObjectPointer ()
    {
      if( object_ != 0 )
        object_->removeReference();
    }

    /** \brief assign another pointer to this one.
     */
    inline ThisType &operator= ( const ThisType &other )
    {
      // Note that it is safe to remove the reference first. If other holds
      // a reference to the same object, the reference counter cannot reach
      // zero.
      if( object_ != 0 )
        object_->removeReference();
      object_ = other.object_;
      if( object_ != 0 )
        object_.addReference();
    }

    /** \brief dereference the ObjectPointer
     *
     *  \note This method asserts that the pointer is not 0.
     * 
     *  \returns a reference to the object pointed to.
     */
    ObjectType &operator* () const
    {
      assert( object_ != 0 );
      return object_->getObject();
    }
  };



#if 0
  /** \class ObjectReference
   *  \brief models a reference to a reference countable object
   *
   *  ObjectReference tries to behave just like a normal C++ reference to an
   *  object implementing the ReferenceCounterInterface. Internally, however,
   *  it calls the object's addReference method when the reference is created
   *  and the removeReference method when it is deleted.
   */
  template< class ObjectImp >
  class ObjectReference
  {
  public:
    //! type of object this reference points to
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
    inline operator const ObjectType& () const
    {
      return object_;
    }

    inline operator ObjectType& ()
    {
      return object_;
    }

    // somehow necessary for LocalFunctionWrapper (at the moment...)
    const ObjectType &operator* () const
    {
      return object_;
    }

    ObjectType &operator* ()
    {
      return object_;
    }
  };
#endif

}

#endif
