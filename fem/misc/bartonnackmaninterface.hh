#include <dune/common/bartonnackmanifcheck.hh>

#ifndef DUNE_FEM_BARTONNACKMANINTERFACE_HH
#define DUNE_FEM_BARTONNACKMANINTERFACE_HH

#include <dune/common/misc.hh>
#include <dune/common/typetraits.hh>

namespace Dune
{
  template< class Interface, class Implementation >
  class BartonNackmanInterface
  {
  private:
    typedef BartonNackmanInterface< Interface, Implementation > ThisType;

  public:
    inline BartonNackmanInterface ()
    {
      // make sure the interface is derived from BartonNackmanInterface
      typedef CompileTimeChecker< Conversion< Interface, ThisType > :: exists >
        __Interface_Must_Be_Derived_From_BartonNackmanInterface__;
      // make sure the implementation is derived from its interface
      typedef CompileTimeChecker< Conversion< Interface, ThisType > :: exists >
        __Implementation_Must_Be_Derived_From_Its_Interface__;
    }
    
  protected:
    inline const Implementation &asImp () const
    {
      return static_cast< const Implementation & >( *this );
    }

    inline Implementation &asImp ()
    {
      return static_cast< Implementation & >( *this );
    }
  };
  
}

#endif
