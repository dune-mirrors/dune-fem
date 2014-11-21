#include <dune/common/bartonnackmanifcheck.hh>

#ifndef DUNE_FEM_BARTONNACKMANINTERFACE_HH
#define DUNE_FEM_BARTONNACKMANINTERFACE_HH

#include <dune/common/typetraits.hh>

namespace Dune
{

  namespace Fem
  {

    template< class Interface, class Implementation >
    class BartonNackmanInterface
    {
      typedef BartonNackmanInterface< Interface, Implementation > ThisType;

    public:
      BartonNackmanInterface ()
      {
        static_assert( (Conversion< Interface, ThisType >::exists), "Interface must be derived from BartonNackmanInterface." );
        //static_assert( (Conversion< Implementation, Interface >::exists), "Implementation must be derived from its interface." );
      }

    protected:
      static const Implementation &asImp ( const ThisType &other )
      {
        return static_cast< const Implementation & >( other );
      }

      static Implementation &asImp ( ThisType &other )
      {
        return static_cast< Implementation & >( other );
      }

      const Implementation &asImp () const
      {
        return asImp( *this );
      }

      Implementation &asImp ()
      {
        return asImp( *this );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_BARTONNACKMANINTERFACE_HH
