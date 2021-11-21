#ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_DATAHANDLE_HH
#define DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_DATAHANDLE_HH

#include <type_traits>

#include <dune/grid/common/datahandleif.hh>

#include <dune/fem/gridpart/geometrygridpart/entity.hh>

namespace Dune
{

  namespace Fem
  {

    // GeometryGridPartDataHandle
    // ------------

    template< class GridFamily, class WrappedHandle >
    class GeometryGridPartDataHandle
    : public CommDataHandleIF< GeometryGridPartDataHandle< GridFamily, WrappedHandle >, typename WrappedHandle::DataType >
    {
      typedef typename std::remove_const< GridFamily >::type::Traits Traits;
      typedef typename GridFamily::GridFunctionType GridFunctionType;

      template< class HostEntity >
      class EntityProxy;

    public:
      GeometryGridPartDataHandle ( WrappedHandle &handle, const GridFunctionType &gridFunction  )
        : wrappedHandle_( handle ), gridFunction_( gridFunction )
      {}

      bool contains ( int dim, int codim ) const
      {
        return wrappedHandle_.contains( dim, codim );
      }

      bool fixedSize ( int dim, int codim ) const
      {
        return wrappedHandle_.fixedSize( dim, codim );
      }

      template< class HostEntity >
      size_t size ( const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity > proxy( gridFunction_, hostEntity );
        return wrappedHandle_.size( *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        EntityProxy< HostEntity > proxy( gridFunction_, hostEntity );
        wrappedHandle_.gather( buffer, *proxy );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
      {
        EntityProxy< HostEntity > proxy( gridFunction_, hostEntity );
        wrappedHandle_.scatter( buffer, *proxy, size );
      }

    private:
      WrappedHandle &wrappedHandle_;
      const GridFunctionType &gridFunction_;
    };



    template< class GridFamily, class WrappedHandle >
    template< class HostEntity >
    class GeometryGridPartDataHandle< GridFamily, WrappedHandle >::EntityProxy
    {
    public:
      static const int dimension = HostEntity::dimension;
      static const int codimension = HostEntity::codimension;

      typedef typename GridFamily::template Codim< codimension > :: Entity  Entity;

    protected:
      typedef typename Entity::Implementation  EntityImpl;

    public:
      EntityProxy ( const GridFunctionType &gridFunction_, const HostEntity &hostEntity )
      : entity_( EntityImpl( gridFunction_, hostEntity ) )
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

#endif // #ifndef DUNE_FEM_GRIDPART_GEOMETRYGRIDPART_DATAHANDLE_HH
