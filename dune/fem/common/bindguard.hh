#ifndef DUNE_FEM_COMMON_BINDGUARD_HH
#define DUNE_FEM_COMMON_BINDGUARD_HH

#include <memory>
#include <utility>

namespace Dune
{

  namespace Fem
  {

    // BindGuard
    // ---------

    template< class Object >
    struct BindGuard
    {
      template< class... Args >
      BindGuard ( Object &object, Args &&... args )
        : object_( &object )
      {
        object.bind( std::forward< Args >( args )... );
      }

    private:
      struct Deleter
      {
        void operator() ( Object *object ) { object->unbind(); }
      };

      std::unique_ptr< Object, Deleter > object_;
    };



    // bindGuard
    // ---------

    template< class Object, class... Args >
    BindGuard< Object > bindGuard ( Object &object, Args &&... args )
    {
      return BindGuard< Object >( object, std::forward< Args >( args )... );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_BINDGUARD_HH
