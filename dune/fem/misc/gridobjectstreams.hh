#ifndef DUNE_FEM_GRIDOBJECTSTREAMS_HH
#define DUNE_FEM_GRIDOBJECTSTREAMS_HH

#include <dune/common/exceptions.hh>

#include <dune/fem/misc/griddeclaration.hh>

namespace Dune
{

  namespace Fem
  {

    // DummyObjectStream
    // -----------------

    struct DummyObjectStream
    {
      class EOFException {};

      template< class T >
      void read ( T & ) const { DUNE_THROW( NotImplemented, "DummyObjectStream::read not implemented." ); }

      template< class T >
      void readObject ( T & ) { DUNE_THROW( NotImplemented, "DummyObjectStream::readObject not implemented." ); }

      void readObject ( int ) { DUNE_THROW( NotImplemented, "DummyObjectStream::readObject not implemented." ); }
      void readObject ( double ) { DUNE_THROW( NotImplemented, "DummyObjectStream::readObject not implemented." ); }

      template< class T >
      void write ( const T & ) { DUNE_THROW( NotImplemented, "DummyObjectStream::write not implemented." ); }

      template< class T >
      void writeObject ( T & ) { DUNE_THROW( NotImplemented, "DummyObjectStream::writeObject not implemented." ); }

      void writeObject ( int ) { DUNE_THROW( NotImplemented, "DummyObjectStream::writeObject not implemented." ); }
      void writeObject ( double ) { DUNE_THROW( NotImplemented, "DummyObjectStream::writeObject not implemented." ); }
    };



    // GridObjectStreamTraits
    // ----------------------

    template< class Grid >
    struct GridObjectStreamTraits
    {
      typedef DummyObjectStream InStreamType;
      typedef DummyObjectStream OutStreamType;
    };

    template< class Grid >
    struct GridObjectStreamTraits< const Grid >
    {
      typedef typename GridObjectStreamTraits< Grid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< Grid >::OutStreamType OutStreamType;
    };



    // GridObjectStreamTraits for ALUGrid
    // ----------------------------------

#if HAVE_DUNE_ALUGRID
    template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct GridObjectStreamTraits< ALUGrid< dim, dimworld, elType, refineType, Comm > >
    {
      typedef typename ALUGrid< dim, dimworld, elType, refineType, Comm >::ObjectStreamType InStreamType;
      typedef typename ALUGrid< dim, dimworld, elType, refineType, Comm >::ObjectStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_ALUGRID



    // GridObjectStreamTraits for CacheItGrid
    // --------------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct GridObjectStreamTraits< CacheItGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< HostGrid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< HostGrid >::OutStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_METAGRID



    // GridObjectStreamTraits for CartesianGrid
    // ----------------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct GridObjectStreamTraits< CartesianGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< HostGrid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< HostGrid >::OutStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_METAGRID



    // GridObjectStreamTraits for FilteredGrid
    // ---------------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct GridObjectStreamTraits< FilteredGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< HostGrid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< HostGrid >::OutStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_METAGRID



    // GridObjectStreamTraits for GeometryGrid
    // ---------------------------------------

    template< class HostGrid, class CoordFunction, class Allocator >
    struct GridObjectStreamTraits< GeometryGrid< HostGrid, CoordFunction, Allocator > >
    {
      typedef typename GridObjectStreamTraits< HostGrid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< HostGrid >::OutStreamType OutStreamType;
    };



    // GridObjectStreamTraits for IdGrid
    // ---------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct GridObjectStreamTraits< IdGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< HostGrid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< HostGrid >::OutStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_METAGRID



    // GridObjectStreamTraits for ParallelGrid
    // ---------------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct GridObjectStreamTraits< ParallelGrid< HostGrid > >
    {
      typedef typename ParallelGrid< HostGrid >::RankManager::ObjectStream InStreamType;
      typedef typename ParallelGrid< HostGrid >::RankManager::ObjectStream OutStreamType;
    };
#endif // #if HAVE_DUNE_METAGRID



    // GridObjectStreamTraits for SphereGrid
    // -------------------------------------

#if HAVE_DUNE_METAGRID
    template< class HostGrid, class MapToSphere >
    struct GridObjectStreamTraits< SphereGrid< HostGrid, MapToSphere > >
    {
      typedef typename GridObjectStreamTraits< HostGrid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< HostGrid >::OutStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_METAGRID

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDOBJECTSTREAMS_HH
