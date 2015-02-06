#ifndef DUNE_IDGRID_DATAHANDLE_HH
#define DUNE_IDGRID_DATAHANDLE_HH

//- dune-common includes
#include <dune/common/typetraits.hh>

//- dune-grid includes
#include <dune/grid/common/datahandleif.hh>

#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/common/ldbhandleif.hh>
#endif

namespace Dune
{

  // IdDataHandle
  // ------------

  template< class WrappedHandle, class GridFamily >
  class IdDataHandle
    : public CommDataHandleIF< IdDataHandle< WrappedHandle, GridFamily >, typename WrappedHandle::DataType >
  {
  protected:
    typedef IdDataHandle< WrappedHandle, GridFamily > ThisType;

    // type of traits
    typedef typename remove_const< GridFamily >::type::Traits Traits;

    typedef typename Traits :: ExtraData ExtraData ;

    template< int codim >
    struct Codim
    {
      // type of entity
      typedef typename Traits::template Codim< codim >::EntityType EntityType;
    };

  public:
    // type of data to be communicated
    typedef typename WrappedHandle::DataType DataType;

    typedef CommDataHandleIF< ThisType, DataType > DataHandleIFType;

  private:
    // prohibit copying
    IdDataHandle ( const ThisType & );

  public:
    IdDataHandle ( ExtraData data, WrappedHandle &wrappedHandle )
    : wrappedHandle_( wrappedHandle ),
      data_( data )
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
      typedef typename Codim< HostEntity::codimension >::EntityType EntityType;
      const EntityType entity( typename EntityType::Implementation( data(), hostEntity ) );
      return wrappedHandle_.size( entity );
    }

    template< class MessageBuffer, class HostEntity >
    void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
    {
      typedef typename Codim< HostEntity::codimension >::EntityType EntityType;
      const EntityType entity( typename EntityType::Implementation( data(), hostEntity ) );
      wrappedHandle_.gather( buffer, entity );
    }

    template< class MessageBuffer, class HostEntity >
    void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
    {
      typedef typename Codim< HostEntity::codimension >::EntityType EntityType;
      const EntityType entity( typename EntityType::Implementation( data(), hostEntity ) );
      wrappedHandle_.scatter( buffer, entity, size );
    }

    ExtraData data() const { return data_; }

  protected:
    WrappedHandle &wrappedHandle_;
    ExtraData data_;
  };

} // namespace Dune

#endif // #ifndef DUNE_IDGRID_DATAHANDLE_HH
