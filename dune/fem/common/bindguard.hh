#ifndef DUNE_FEM_COMMON_BINDGUARD_HH
#define DUNE_FEM_COMMON_BINDGUARD_HH

#include <memory>
#include <tuple>
#include <utility>

#include <dune/common/typeutilities.hh>

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



    // isBindable
    // ----------

    namespace Impl
    {

      template< class Object, class... Args >
      auto isBindable ( PriorityTag< 1 >, Object &object, Args &&... args )
        -> decltype( std::declval< Object & >().bind( std::declval< Args >()... ), std::true_type() );

      template< class Object, class... Args >
      auto isBindable ( PriorityTag< 0 >, Object &object, Args &&... args )
        -> std::false_type;

    } // namespace Impl

    template< class Object, class... Args >
    struct isBindable
      : public decltype( Impl::isBindable( PriorityTag< 42 >(), std::declval< Object & >(), std::declval< Args >()... ) )
    {};



    // bindGuard
    // ---------

    template< class Object, class... Args >
    inline static auto bindGuard ( Object &object, Args &&... args )
      -> std::enable_if_t< isBindable< Object, Args... >::value, BindGuard< Object > >
    {
      return BindGuard< Object >( object, std::forward< Args >( args )... );
    }

    template< std::size_t... i, class Objects, class... Args >
    inline static auto bindGuard ( std::index_sequence< i... >, Objects objects, Args &&... args )
      -> std::tuple< decltype( bindGuard( std::declval< std::tuple_element_t< i, Objects > >(), std::declval< Args >()... ) )... >
    {
      // note: do not forward args, as there the tuple might contain multiple entries
      return std::make_tuple( bindGuard( std::get< i >( objects ), args... )... );
    }

    template< class... Object, class... Args >
    inline static auto bindGuard ( std::tuple< Object &... > objects, Args &&... args )
      -> decltype( bindGuard( std::index_sequence_for< Object... >(), objects, std::forward< Args >( args )... ) )
    {
      return bindGuard( std::index_sequence_for< Object... >(), objects, std::forward< Args >( args )... );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_BINDGUARD_HH
