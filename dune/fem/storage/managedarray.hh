#ifndef DUNE_FEM_HPDG_UTILITY_MANAGEDARRAY_HH
#define DUNE_FEM_HPDG_UTILITY_MANAGEDARRAY_HH

#include <cstddef>

#include <algorithm>
#include <tuple>
#include <vector>

#include <dune/fem/space/common/dofmanager.hh>
#include <dune/fem/space/common/restrictprolonginterface.hh>
#include <dune/fem/space/mapper/codimensionmapper.hh>

namespace Dune
{

  namespace Fem
  {

    namespace hpDG
    {

      // ManagedArray
      // ------------

      /** \brief An associative container that may be restricted and prolonged during local mesh adaptation
       *
       *  \tparam GridPart  a Dune::Fem::GridPart
       *
       *  \ingroup DiscreteFunctionSpace_API
       */
      template< class GridPart, class T = int >
      class ManagedArray
      {
        using container = std::vector< T >;

      public:
        /** \brief value type */
        using value_type = typename container::value_type;

        /** \brief size type */
        using size_type = typename container::size_type;

        /** \brief iterator type */
        using iterator = typename container::iterator;
        /** \brief const iterator type */
        using const_iterator = typename container::const_iterator;

        /** \brief reference type */
        using reference = typename container::reference;
        /** \brief const reference type */
        using const_reference = typename container::const_reference;

        /** \name Construction
         *  \{
         */

        ManagedArray ( GridPart &gridPart, const value_type& p )
          : mapper_( gridPart )
        {
          std::tie( dofStorage_, p_ )
            = Dune::Fem::allocateManagedDofStorage( gridPart.grid(), mapper_, "ManagedArray", p_ );
          std::fill( begin(), end(), p );
        }

        ManagedArray ( const ManagedArray & ) = delete;

        ManagedArray &operator= ( const ManagedArray & ) = delete;

        /** \} */

        /** \name Iterators
         *  \{
         */

        /** \brief return begin iterator */
        iterator begin () noexcept { return p_->begin(); }

        /** \brief return begin iterator */
        const_iterator begin () const noexcept { return p_->begin(); }

        /** \brief return end iterator */
        iterator end () noexcept { return p_->end(); }

        /** \brief return end iterator */
        const_iterator end () const noexcept { return p_->end(); }

        /** \brief return begin iterator */
        const_iterator cbegin () const noexcept { return p_->cbegin(); }

        /** \brief return end iterator */
        const_iterator cend () const noexcept { return p_->cend(); }

        /** \} */

        /** \name Capacity
         *  \{
         */

        /** \brief return number of elements */
        size_type size () const noexcept { return p_->size(); }

        /** \brief resize array */
        void resize ( size_type size ) { p_->resize( size ); }

        /** \brief reserve memory */
        void reserve ( size_type size ) { p_->reserve( size ); }

        /** \} */

        /** \name Element access
         *  \{
         */

        /** \brief return \f$i\f$-th element */
        reference operator[] ( size_type i ) { return (*p_)[ i ]; }

        /** \brief return data for given element */
        reference operator[] ( const typename GridPart::template Codim< 0 >::EntityType &entity )
        {
          return (*p_)[ index( entity ) ];
        }

        /** \brief return \f$i\f$-th element */
        const_reference operator[] ( size_type i ) const { return (*p_)[ i ]; }

        /** \brief return data for given element */
        const_reference operator[] ( const typename GridPart::template Codim< 0 >::EntityType &entity ) const
        {
          return (*p_)[ index( entity ) ];
        }

        /** \} */

      private:
        friend class Dune::Fem::RestrictProlongDefault< ManagedArray< GridPart > >;

        void enableDofCompression ()
        {
          if( dofStorage_ )
            dofStorage_->enableDofCompression();
        }

        std::size_t index ( const typename GridPart::template Codim< 0 >::EntityType &entity ) const
        {
          std::size_t index;
          mapper_.mapEachEntityDof( entity, [&index]( std::size_t, std::size_t i ){ index = i; } );
          return index;
        }

        Dune::Fem::CodimensionMapper< GridPart, 0 > mapper_;
        Dune::Fem::DofStorageInterface *dofStorage_ = nullptr;
        std::vector< int > *p_ = nullptr;
      };

    } // namespace hpDG



#ifndef DOXYGEN

    // RestrictProlongDefault
    // ----------------------

    template< class GridPart >
    class RestrictProlongDefault< hpDG::ManagedArray< GridPart > >
      : public RestrictProlongInterface< RestrictProlongTraits< RestrictProlongDefault< hpDG::ManagedArray< GridPart > >, typename GridPart::ctype > >
    {
      using ManagedArrayType = hpDG::ManagedArray< GridPart >;

    public:
      explicit RestrictProlongDefault ( ManagedArrayType &p )
        : p_( p )
      {
        p_.enableDofCompression();
      }

      void setFatherChildWeight ( typename GridPart::ctype weight ) const {}

      template< class Entity >
      void restrictLocal ( const Entity &father, const Entity &child, bool initialize ) const
      {
        if( p_.index( father ) == p_.index( child ) )
          return;

        if( initialize )
          p_[ father ] = p_[ child ] + 1;
        else
          p_[ father ] = std::max( p_[ father ], p_[ child ] + 1 );
      }

      template< class Entity, class LocalGeometry >
      void restrictLocal ( const Entity &father, const Entity &child, const LocalGeometry &, bool initialize ) const
      {
        restrictLocal( father, child, initialize );
      }

      template< class Entity >
      void prolongLocal ( const Entity &father, const Entity &child, bool initialize ) const
      {
        p_[ child ] = p_[ father ];
      }

      template< class Entity, class LocalGeometry >
      void prolongLocal ( const Entity &father, const Entity &child, const LocalGeometry &, bool initialize ) const
      {
        prolongLocal( father, child, initialize );
      }

      template< class Communicator >
      void addToList ( Communicator & ) {}

      template< class LoadBalancer >
      void addToLoadBalancer ( LoadBalancer & ) {}

    private:
      ManagedArrayType &p_;
    };

#endif // #ifndef DOXYGEN

    using hpDG::ManagedArray;

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_HPDG_UTILITY_MANAGEDARRAY_HH
