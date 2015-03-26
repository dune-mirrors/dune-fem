#ifndef DUNE_FEM_GRIDPART_COMMON_DEFAULTENTITYPOINTER_HH
#define DUNE_FEM_GRIDPART_COMMON_DEFAULTENTITYPOINTER_HH

#include <utility>

namespace Dune
{

  namespace Fem
  {

    // DefaultEntityPointer
    // --------------------

    template< class E >
    class DefaultEntityPointer
    {
    public:
      /** \brief entity type */
      typedef E Entity;

      /** \brief codimension */
      static const int codimension = Entity::codimension;

      /** \name Construction
       *  \{
       */

      DefaultEntityPointer () = default;

      template< class... Args >
      DefaultEntityPointer ( Args &&... args )
        : entity_( std::forward< Args >( args )... )
      {}

#ifndef DOXYGEN

      // forward to default copy constructor
      DefaultEntityPointer ( DefaultEntityPointer &other )
        : DefaultEntityPointer( static_cast< const DefaultEntityPointer & >( other ) )
      {}

#endif // #ifndef DOXYGEN

      /** \} */

#ifndef DOXYGEN

      DefaultEntityPointer ( const DefaultEntityPointer & ) = default;

      DefaultEntityPointer ( DefaultEntityPointer && ) = default;

      DefaultEntityPointer &operator= ( const DefaultEntityPointer & ) = default;

      DefaultEntityPointer &operator= ( DefaultEntityPointer && ) = default;

#endif // #ifndef DOXYGEN

      /** \name Public member methods
       *  \{
       */

      Entity dereference () const
      {
        return entity_;
      }

      int level () const
      {
        return entity_.level();
      }

      /** \} */

    private:
      Entity entity_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_COMMON_DEFAULTENTITYPOINTER_HH
