#ifndef DUNE_FEM_ENGINE_HH
#define DUNE_FEM_ENGINE_HH

#ifndef DUNE_FEM_COMPATIBILITY
#define DUNE_FEM_COMPATIBILITY 1
#endif

namespace Dune
{

  template< class Impl, class User >
  class EngineWrapper
  {
  protected:
    inline const Impl &asImp () const
    {
      const User &user = static_cast< const User & >( *this );
      return user.asImp();
    }

    inline Impl &asImp ()
    {
      User &user = static_cast< User & >( *this );
      return user.asImp();
    }
  };



  template< class Impl >
  class EngineDefault
  {
  protected:
    inline const Impl &asImp () const
    {
      return static_cast< const Impl & >( *this );
    }

    inline Impl &asImp ()
    {
      return static_cast< Impl & >( *this );
    }
  };
  
}

#endif
