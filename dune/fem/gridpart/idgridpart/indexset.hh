#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH

#include <type_traits>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/fem/gridpart/common/indexset.hh>
#include <dune/fem/gridpart/common/persistentindexset.hh>

#include <dune/fem/io/streams/streams.hh>

namespace Dune
{

  namespace Fem
  {

    // Internal forward declaration
    // ----------------------------

    template< class GridFamily >
    class IdIndexSet;



    namespace __IdIndexSet
    {

      // IndexSet
      // --------

      template< class GridFamily >
      class IndexSet
      {
      protected:
        typedef typename std::remove_const< GridFamily >::type::Traits Traits;

      public:
        typedef typename Traits::HostGridPartType::IndexSetType HostIndexSetType;

        static const int dimension = HostIndexSetType::dimension;

        template< int codim >
        struct Codim
        {
          typedef typename std::remove_const< GridFamily >::type::Traits::template Codim< codim >::Entity Entity;
        };

        typedef typename HostIndexSetType::IndexType IndexType;

        typedef typename HostIndexSetType::Types Types;

        explicit IndexSet ( const HostIndexSetType &hostIndexSet )
          : hostIndexSet_( &hostIndexSet )
        {}

        Types types ( int codim ) const
        {
          return hostIndexSet().types( codim );
        }

        const std::vector< GeometryType > &geomTypes ( int codim ) const
        {
          return hostIndexSet().geomTypes( codim );
        }

        template< class Entity >
        bool contains ( const Entity &entity ) const
        {
          return hostIndexSet().contains( entity.impl().hostEntity() );
        }

        IndexType size ( GeometryType type ) const
        {
          return hostIndexSet().size( type );
        }

        IndexType size ( int codim ) const
        {
          return hostIndexSet().size( codim );
        }

        template< class Entity >
        IndexType index ( const Entity &entity ) const
        {
          return index< Entity::codimension >( entity );
        }

        template< int codim >
        IndexType index ( const typename Codim< codim >::EntityType &entity ) const
        {
          return hostIndexSet().template index< codim >( entity.impl().hostEntity() );
        }

        template< class Entity >
        IndexType subIndex ( const Entity &entity, int i, unsigned int cd ) const
        {
          return subIndex< Entity::codimension >( entity, i, cd );
        }

        template< int codim >
        IndexType subIndex ( const typename Codim< codim >::Entity &entity, int i, unsigned int cd ) const
        {
          return hostIndexSet().template subIndex< codim >( entity.impl().hostEntity(), i, cd );
        }

        const HostIndexSetType &hostIndexSet () const
        {
          assert( hostIndexSet_ );
          return *hostIndexSet_;
        }

      private:
        const HostIndexSetType &hostIndexSet_;
      };



      // ConsecutiveIndexSet
      // -------------------

      template< class GridFamily >
      class ConsecutiveIndexSet
        : public IndexSet< GridFamily >
      {
        typedef IndexSet< GridFamily > BaseType;

      public:
        typedef typename BaseType::HostIndexSetType HostIndexSetType;

        using BaseType::BaseType;
        using BaseType::hostIndexSet;

        bool consecutive () const { return hostIndexSet().consecutive(); }

        void resize () { hostIndexSet().resize(); }

        bool compress () { return hostIndexSet().compress(); }

        template< class Entity >
        void insertEntity ( const Entity &entity )
        {
          insertEntity< Entity::codimension >( entity );
        }

        template< int codim >
        void insertEntity ( const typename BaseType::template Codim< codim >::EntityType &entity )
        {
          hostIndexSet().template insertEntity< codim >( entity.impl().hostEntity() );
        }

        template< class Entity >
        void removeEntity ( const Entity &entity )
        {
          removeEntity< Entity::codimension >( entity );
        }

        template< int codim >
        void removeEntity ( const typename BaseType::template Codim< codim >::EntityType &entity )
        {
          hostIndexSet().template removeEntity< codim >( entity.impl().hostEntity() );
        }

        void backup () const { hostIndexSet().backup(); }

        void restore () { hostIndexSet().restore(); }

        template< class T >
        void write ( OutStreamInterface< T > &stream )
        {
          hostIndexSet().write( stream );
        }

        template< class T >
        void read ( InStreamInterface< T > &stream )
        {
          hostIndexSet().read( stream );
        }

      protected:
        HostIndexSetType &hostIndexSet ()
        {
          return const_cast< HostIndexSetType & >( BaseType::hostIndexSet() );
        }
      };



      // AdaptiveIndexSet
      // ----------------

      template< class GridFamily >
      class AdaptiveIndexSet
        : public ConsecutiveIndexSet< GridFamily >
      {
        typedef ConsecutiveIndexSet< GridFamily > BaseType;

      public:
        using BaseType::BaseType;

        int numberOfHoles ( GeometryType type ) const
        {
          return this->hostIndexSet().numberOfHoles( type );
        }
        
        int oldIndex ( int hole, GeometryType type ) const
        {
          return this->hostIndexSet().oldIndex( hole, type );
        }
        
        int newIndex ( int hole, GeometryType type ) const
        {
          return this->hostIndexSet().newIndex( hole, type );
        }
      };



      // Implementation
      // --------------

      template< class GridFamily,
                class HostIndexSet = typename std::remove_const< GridFamily >::type::Traits::HostGridPartType::IndexSetType,
                bool consecutive = Capabilities::isConsecutiveIndexSet< HostIndexSet >::v,
                bool adaptive = Capabilities::isAdaptiveIndexSet< HostIndexSet >::v >
      class Implementation
      {
        typedef typename std::conditional< adaptive,
            AdaptiveIndexSet< GridFamily >,
            typename std::conditional< consecutive,
                ConsecutiveIndexSet< GridFamily >,
                IndexSet< GridFamily >
              >::type
          >::type Type;
      };

    } // namespace __IdIndexSet



    // IdIndexSet
    // ----------

    template< class GridFamily >
    class IdIndexSet
      : public __IdIndexSet::Implementation< GridFamily >::Type
    {
      typedef typename __IdIndexSet::Implementation< GridFamily >::Type BaseType;

    public:
      using BaseType::BaseType;
    };



    namespace Capabilities
    {

      template< class GridFamily >
      struct isConsecutiveIndexSet< IdIndexSet< GridFamily > >
        : public isConsecutiveIndexSet< typename IdIndexSet< GridFamily >::HostIndexSetType >
      {};

      template< class GridFamily >
      struct isAdaptiveIndexSet< IdIndexSet< GridFamily > >
        : public isAdaptiveIndexSet< typename IdIndexSet< GridFamily >::HostIndexSetType >
      {};

      template< class GridFamily >
      struct isPersistentIndexSet< IdIndexSet< GridFamily > >
      {
      private:
        typedef IdIndexSet< GridFamily > IndexSetType;
        typedef typename IndexSetType::HostIndexSetType HostIndexSetType;

      public:
        static const bool v = isPersistentIndexSet< HostIndexSetType >::v;

        static constexpr PersistentIndexSetInterface *map ( IndexSetType &indexSet ) noexcept
        {
          return isPersistentIndexSet< HostIndexSetType >::map( indexSet.hostIndexSet() );
        }
      };

    } // namespace Capabilities

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH
