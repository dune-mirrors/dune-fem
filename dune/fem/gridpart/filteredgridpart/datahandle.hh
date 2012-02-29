#ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_DATAHANDLE_HH
#define DUNE_FEM_GRIDPART_FILTEREDGRIDPART_DATAHANDLE_HH


//- dune-common includes
#include <dune/common/typetraits.hh>

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
    {
      typedef GridPart GridPartType;
      typedef typename remove_const< GridPartType >::type::Traits Traits;

    public:
      FilteredGridPartDataHandle ( WrappedHandle &dataHandle, const GridPart &gridPart )
      : gridPart_( gridPart ),
        wrappedHandle_( dataHandle )
      { } 

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

  }  // end namespace Fem

}  // end namespace Dune

#endif // #ifndef DUNE_FEM_GRIDPART_FILTEREDGRIDPART_DATAHANDLE_HH
