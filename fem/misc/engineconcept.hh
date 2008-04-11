#ifndef DUNE_FEM_ENGINE_HH
#define DUNE_FEM_ENGINE_HH

namespace Dune
{

  /** \addtogroup EngineConcept
   *
   *  The engine concept is a way to define static interfaces. In an engine
   *  concept, the interface is a wrapper around the actual implementation
   *  (called engine in this context). In contrast to Barton-Nackman
   *  interfaces, the user only deals with the interface class.
   *
   *  Still, the classical engine concept lacks two things:
   *  - It is not possible to provide a default implementation with the same
   *    power as for Barton-Nackman interfaces.
   *  - The implementation actually resides within the interface class.
   *  .
   *  To overcome these difficulties, we combine the power of both interface
   *  types. The actual interface uses the Barton-Nackman trick to obtain the
   *  implementation for a socalled \b User class. This class has two purposes:
   *  - It knows how to construct the implementation, so it can declare the
   *    necessary constructors.
   *  - It knows, where the implementation resides in memory, so it can return
   *    a reference. This is done via the method asImp.
   *  .
   *  To obtain default implementations, we use the Barton-Nackman trick again
   *  (in the usual manner). The advantage over the standard Barton-Nackman
   *  trick is that a default implemented method never calls itself. So, if an
   *  implementation does not implement a method that is not default
   *  implemented, the compiler will throw an error.
   *
   *  To use the engine concept, just do the following:
   *  - derive the interface class from EngineWrapper and tell the compiler
   *    that you are using EngineWrapper::asImp
   *  - for any interface method, use asImp to call the method on the
   *    implementation
   *  - derive the default implementation from EngineDefault and tell the
   *    compiler that you are using EngineDefault::asImp
   *  - to call another function in the default implementation, use asImp, so
   *    that the implementation has the possibility to intercept this call
   *  - for each implementation write an extra user object, which is derived
   *    from the interface and provides the asImp-method
   *  - mark EngineWrapper a friend in any user class
   *  .
   *
   *  \warning The Barton-Nackman trick has the big disadvantage of allowing
   *           for loops, if the interface method is not implemented. In our
   *           engine concept, this can only happen with the asImp method on
   *           the interface. So, when writing user classees, be careful to
   *           implement \b both asImp methods.
   */

  /** \class EngineWrapper
   *  \ingroup EngineConcept
   *
   *  \brief base class for interfaces (engine wrappers)
   *
   *  This class uses the Barton-Nackman trick to privide the asImp method to
   *  the interface class (which is derived from this class).
   *
   *  \note Do not forget to mark this class a friend in any user class.
   *
   *  \warning The Barton-Nackman trick has the big disadvantage of allowing
   *           for loops, if the interface method is not implemented. So, when
   *           writing user classees, be careful to implement \b both asImp
   *           methods.
   *
   *  For more details, we refer to the \ref EngineConcept "Engine Concept".
   */
  template< class Impl, class User >
  class EngineWrapper
  {
  protected:
    EngineWrapper ()
    {
    }
    
  private:
    // Prohibit automatic copying of the interface
    EngineWrapper ( const EngineWrapper & );
    
  protected:
    /** \brief obtain the implementation from the user class (const version)
     *
     *  Ths method uses the Barton-Nackman trick, to obtain a reference to the
     *  implementation from the user class.
     *
     *  \returns a const reference to the implementation
     */
    inline const Impl &asImp () const
    {
      const User &user = static_cast< const User & >( *this );
      return user.asImp();
    }

    /** \brief obtain the implementation from the user class (non-const version)
     *
     *  Ths method uses the Barton-Nackman trick, to obtain a reference to the
     *  implementation from the user class.
     *
     *  \returns a reference to the implementation
     */
    inline Impl &asImp ()
    {
      User &user = static_cast< User & >( *this );
      return user.asImp();
    }
  };



  /** \class EngineDefault
   *  \ingroup EngineConcept
   *
   *  \brief base class of default implementations
   *
   *  This class uses the Barton-Nackman trick to privide the asImp method to
   *  the default implementation (which is derived from this class).
   *
   *  For more details, we refer to the \ref EngineConcept "Engine Concept".
   */
  template< class Impl >
  class EngineDefault
  {
  protected:
    /** \brief obtain the implementation (const version)
     *
     *  Ths method uses the Barton-Nackman trick, to obtain a reference to the
     *  implementation.
     *
     *  \returns a const reference to the implementation
     */
    inline const Impl &asImp () const
    {
      return static_cast< const Impl & >( *this );
    }

    /** \brief obtain the implementation (non-const version)
     *
     *  Ths method uses the Barton-Nackman trick, to obtain a reference to the
     *  implementation.
     *
     *  \returns a reference to the implementation
     */
    inline Impl &asImp ()
    {
      return static_cast< Impl & >( *this );
    }
  };
  
}

#endif
