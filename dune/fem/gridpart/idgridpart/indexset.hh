#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH

#include <type_traits>
#include <vector>

#include <dune/geometry/type.hh>

#include <dune/fem/gridpart/common/capabilities.hh>
#include <dune/fem/gridpart/common/persistentindexset.hh>

namespace Dune
{

  namespace Fem
  {

    // IdIndexSet
    // ----------

    template< class GridFamily,
              bool isAdaptive = GridPartCapabilities::hasAdaptiveIndexSet< typename std::remove_const< GridFamily >::type::Traits::HostGridPartType >::v >
    class IdIndexSet;



    // IdIndexSetBase
    // --------------

    template< class GridFamily >
    class IdIndexSetBase
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;
      typedef typename Traits::HostGridPartType HostGridPartType;

    protected:
      typedef typename HostGridPartType::IndexSetType HostIndexSetType;

    public:
      static const int dimension = HostIndexSetType::dimension;

      template< int codim >
      struct Codim
      {
        typedef typename std::remove_const< GridFamily >::type::Traits::template Codim< codim >::Entity Entity;
      };

      typedef typename HostIndexSetType::IndexType IndexType;

      typedef typename HostIndexSetType::Types Types;

      explicit IdIndexSetBase ( const HostIndexSetType &hostIndexSet )
        : hostIndexSet_( const_cast< HostIndexSetType * > (&hostIndexSet) )
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

      IndexType size ( const GeometryType &type ) const
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
        return hostIndexSet().subIndex( entity.impl().hostEntity(), i, cd );
      }

    protected:
      const HostIndexSetType &hostIndexSet () const
      {
        assert( hostIndexSet_ );
        return *hostIndexSet_;
      }

      HostIndexSetType &hostIndexSet ()
      {
        assert( hostIndexSet_ );
        return *hostIndexSet_;
      }

    private:
      HostIndexSetType *hostIndexSet_;
    };



    // IdIndexSet
    // ----------

    template< class GridFamily >
    class IdIndexSet< GridFamily, false >
      : public IdIndexSetBase< GridFamily >
    {
      typedef IdIndexSetBase< GridFamily > BaseType;

      template< class > friend class isPersistentIndexSet;

    public:
      using BaseType::BaseType;
    };

    template< class GridFamily >
    class IdIndexSet< GridFamily, true >
      : public IdIndexSetBase< GridFamily >
    {
      typedef IdIndexSetBase< GridFamily > BaseType;

      template< class > friend class isPersistentIndexSet;

    protected:
      using BaseType::hostIndexSet;

    public:
      using BaseType::BaseType;

      bool consecutive () const
      {
        return hostIndexSet().consecutive();
      }

      void resize () { hostIndexSet.resize(); }

      void compress ()
      {
        // assert( consecutive() );
        hostIndexSet().compress();
      }

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

      int numberOfHoles ( GeometryType type ) const
      {
        return hostIndexSet().numberOfHoles( type );
      }

      int oldIndex ( int hole, GeometryType type ) const
      {
        return hostIndexSet().oldIndex( hole, type );
      }

      int newIndex ( int hole, GeometryType type ) const
      {
        return hostIndexSet().newIndex( hole, type );
      }
    };



    // Tempate specialization of isPersistentIndexSet
    // ----------------------------------------------

    template< class GridFamily, bool isAdaptive >
    struct isPersistentIndexSet< IdIndexSet< GridFamily, isAdaptive > >
    {
    private:
      typedef IdIndexSet< GridFamily, isAdaptive > IndexSetType;
      typedef typename IndexSetType::HostIndexSetType HostIndexSetType;

    public:
      static const bool v = isPersistentIndexSet< HostIndexSetType >::v;

      static constexpr PersistentIndexSetInterface *map ( IndexSetType &indexSet ) noexcept
      {
        return isPersistentIndexSet< HostIndexSetType >::map( indexSet.hostIndexSet() );
      }
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH
