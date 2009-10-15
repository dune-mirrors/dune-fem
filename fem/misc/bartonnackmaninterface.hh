#include <dune/common/bartonnackmanifcheck.hh>

#ifndef DUNE_FEM_BARTONNACKMANINTERFACE_HH
#define DUNE_FEM_BARTONNACKMANINTERFACE_HH

#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

namespace Dune
{
  template< class Interface, class Implementation >
  class BartonNackmanInterface
  {
    typedef BartonNackmanInterface< Interface, Implementation > ThisType;

  public:
    BartonNackmanInterface ()
    {
      dune_static_assert( (Conversion< Interface, ThisType >::exists), "Interface must be derived from BartonNackmanInterface." );
      dune_static_assert( (Conversion< Interface, ThisType >::exists), "Implementation must be derived from its interface." );
    }
    
  protected:
    const Implementation &asImp () const
    {
      return static_cast< const Implementation & >( *this );
    }

    Implementation &asImp ()
    {
      return static_cast< Implementation & >( *this );
    }
  };
  
}

#endif // #ifndef DUNE_FEM_BARTONNACKMANINTERFACE_HH
