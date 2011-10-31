#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH

#include <vector>

#include <dune/common/typetraits.hh>

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/indexidset.hh>

namespace Dune
{

  namespace Fem
  {

    // IdIndexSet
    // ----------

    template< class HostIndexSet >
    class IdIndexSet
    {
      typedef IdIndexSet< HostIndexSet > This;

    public:
      typedef typename HostIndexSet::IndexType IndexType;

      IdIndexSet ( const HostIndexSet &hostIndexSet )
      : hostIndexSet_( &hostIndexSet )
      {}

      int size ( const GeometryType &type ) const { return hostIndexSet().size( type ); }
      int size ( const int codim ) const { return hostIndexSet().size( codim ); }

      template< class Entity >
      int index ( const Entity &entity ) const
      {
        return hostIndexSet().index( entity.impl().hostEntity() );
      }

      template< class Entity >
      int subIndex ( const Entity &entity, const int local, const unsigned int codim ) const
      {
        return hostIndexSet().subIndex( entity.impl().hostEntity(), local, codim );
      }

      const std::vector< GeometryType > &geomTypes ( const int codim ) const
      {
        return hostIndexSet().geomTypes( codim );
      }

      template< class Entity >
      bool contains ( const Entity &entity ) const
      {
        static const int cc = Entity::codimension;
        return hostIndexSet().contains( entity.impl().hostEntity() );
      }

      bool adaptive () const { return hostIndexSet().adaptive(); }
      bool persistent () const { return hostIndexSet().persistent(); }

      int numberOfHoles ( const int codim ) const { return hostIndexSet().numberOfHoles( codim ); }

      int oldIndex ( const int hole, const int codim ) const
      {
        return hostIndexSet().oldIndex( hole, codim );
      }

      int newIndex ( const int hole, const int codim ) const
      {
        return hostIndexSet().newIndex( hole, codim );
      }

    private:
      const HostIndexSet &hostIndexSet () const
      {
        assert( hostIndexSet_ );
        return *hostIndexSet_;
      }

      const HostIndexSet *hostIndexSet_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_INDEXSET_HH
