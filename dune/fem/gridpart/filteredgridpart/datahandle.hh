#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_DATAHANDLE_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_DATAHANDLE_HH

// C++
#include <type_traits>

//- dune-grid includes
#include <dune/grid/common/datahandleif.hh>


namespace Dune
{

  namespace Fem
  {

    // Forward declaration
    // -------------------

    template< class HostGridPartImp, class FilterImp, bool useFilteredIndexSet >
    class FilteredGridPart;



    // FilteredGridPartDataHandle
    // --------------------------

    template< class WrappedHandle, class GridPart >
    class FilteredGridPartDataHandle
    : public CommDataHandleIF< FilteredGridPartDataHandle< WrappedHandle, GridPart >, typename WrappedHandle::DataType >
    {
      typedef CommDataHandleIF< FilteredGridPartDataHandle< WrappedHandle, GridPart >,
              typename WrappedHandle::DataType > BaseType;
      typedef GridPart GridPartType;
      typedef typename std::remove_const< GridPartType >::type::Traits Traits;

    public:
      FilteredGridPartDataHandle ( WrappedHandle &dataHandle, const GridPart &gridPart )
      : gridPart_( gridPart ),
        wrappedHandle_( dataHandle )
      { }

      bool contains ( int dim, int codim ) const
      {
        return wrappedHandle_.contains( dim, codim );
      }

      bool fixedSize ( int dim, int codim ) const
      {
        return false;
      }

      template< class HostEntity >
      size_t size ( const HostEntity &hostEntity ) const
      {
        if( gridPart().contains( hostEntity ) )
          return wrappedHandle_.size( hostEntity );
        else
          return 0;
      }

      template< class MessageBuffer, class HostEntity >
      void gather ( MessageBuffer &buffer, const HostEntity &hostEntity ) const
      {
        if( gridPart().contains( hostEntity ) )
          wrappedHandle_.gather( buffer, hostEntity );
      }

      template< class MessageBuffer, class HostEntity >
      void scatter ( MessageBuffer &buffer, const HostEntity &hostEntity, size_t size )
      {
        if( gridPart().contains( hostEntity ) )
          wrappedHandle_.scatter( buffer, hostEntity, size );
        else
        {
          typename BaseType::DataType tmp;
          for (size_t i=0;i<size;++i)
            buffer.read(tmp);
        }
      }

    protected:
      const GridPart &gridPart () const
      {
        return gridPart_;
      }

    private:
      const GridPart &gridPart_;
      WrappedHandle &wrappedHandle_;
    };

  }  // namespace Fem

}  // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_DATAHANDLE_HH
