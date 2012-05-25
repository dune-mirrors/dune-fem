#ifndef DUNE_FEM_GRIDPART_IDGRIDPART_DATAHANDLE_HH
#define DUNE_FEM_GRIDPART_IDGRIDPART_DATAHANDLE_HH

#include <dune/common/typetraits.hh>

#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/gridpart/idgridpart/entity.hh>

namespace Dune
{

  namespace Fem
  {

    // IdDataHandle
    // ------------

    template< class GridFamily, class WrappedHandle >
    class IdDataHandle
    : public CommDataHandleIF< IdDataHandle< GridFamily, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef typename remove_const< GridFamily >::type::Traits Traits;

      template< class HostEntity >
      class EntityProxy;

    public:
      IdDataHandle ( WrappedHandle &handle )
      : wrappedHandle_( handle )
      {}

      bool contains ( int dim, int codim ) const
      {
        return wrappedHandle_.contains( dim, codim );
      }

      bool fixedsize ( int dim, int codim ) const
      {
        return wrappedHandle_.fixedsize( dim, codim );
      }

      template< class HostEntity >
      size_t size ( const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity > proxy( hostEntity );
        return wrappedHandle_.size( *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity > proxy( hostEntity );
        wrappedHandle_.gather( buffer, *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
      {
        EntityProxy< HostEntity > proxy( hostEntity );
        wrappedHandle_.scatter( buffer, *proxy, size );
      }

    private:
      WrappedHandle &wrappedHandle_;
    };



    template< class GridFamily, class WrappedHandle >
    template< class HostEntity >
    struct IdDataHandle< GridFamily, WrappedHandle >::EntityProxy
    {
      static const int dimension = HostEntity::dimension;
      static const int codimension = HostEntity::codimension;
      
      typedef Dune::Entity< codimension, dimension, const GridFamily, IdEntity > Entity;

    private:
      typedef IdEntity< codimension, dimension, const GridFamily > EntityImpl;

    public:
      EntityProxy ( const HostEntity &hostEntity )
      : entity_( EntityImpl( hostEntity ) )
      {}

      const Entity &operator* () const
      {
        return entity_;
      }

    private:
      Entity entity_;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_IDGRIDPART_DATAHANDLE_HH
