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

#if HAVE_ALUGRID
    template< int dim, int dimworld >
    struct GridObjectStreamTraits< ALUConformGrid< dim, dimworld> >
    {
      typedef typename ALUConformGrid< dim, dimworld >::ObjectStreamType InStreamType;
      typedef typename ALUConformGrid< dim, dimworld >::ObjectStreamType OutStreamType;
    };

    template< int dim, int dimworld >
    struct GridObjectStreamTraits< ALUCubeGrid< dim, dimworld> >
    {
      typedef typename ALUCubeGrid< dim, dimworld >::ObjectStreamType InStreamType;
      typedef typename ALUCubeGrid< dim, dimworld >::ObjectStreamType OutStreamType;
    };

    template< int dim, int dimworld >
    struct GridObjectStreamTraits< ALUSimplexGrid< dim, dimworld> >
    {
      typedef typename ALUSimplexGrid< dim, dimworld >::ObjectStreamType InStreamType;
      typedef typename ALUSimplexGrid< dim, dimworld >::ObjectStreamType OutStreamType;
    };
#endif // #if HAVE_ALUGRID

#if HAVE_ALUGRID || HAVE_DUNE_ALUGRID
    template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
    struct GridObjectStreamTraits< ALUGrid< dim, dimworld, elType, refineType, Comm > >
    {
      typedef typename ALUGrid< dim, dimworld, elType, refineType, Comm >::ObjectStreamType InStreamType;
      typedef typename ALUGrid< dim, dimworld, elType, refineType, Comm >::ObjectStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_ALUGRID


#if HAVE_DUNE_METAGRID
    template< class HostGrid >
    struct GridObjectStreamTraits< CacheItGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< CacheItGrid< HostGrid > >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< CacheItGrid< HostGrid > >::OutStreamType OutStreamType;
    };

    template< class HostGrid >
    struct GridObjectStreamTraits< CartesianGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< CartesianGrid< HostGrid > >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< CartesianGrid< HostGrid > >::OutStreamType OutStreamType;
    };

    template< class HostGrid >
    struct GridObjectStreamTraits< FilteredGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< FilteredGrid< HostGrid > >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< FilteredGrid< HostGrid > >::OutStreamType OutStreamType;
    };

    template< class HostGrid >
    struct GridObjectStreamTraits< IdGrid< HostGrid > >
    {
      typedef typename GridObjectStreamTraits< IdGrid< HostGrid > >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< IdGrid< HostGrid > >::OutStreamType OutStreamType;
    };

    template< class HostGrid >
    struct GridObjectStreamTraits< ParallelGrid< HostGrid > >
    {
      typedef typename ParallelGrid< HostGrid >::RankManager::ObjectStream InStreamType;
      typedef typename ParallelGrid< HostGrid >::RankManager::ObjectStream OutStreamType;
    };

    template< class HostGrid, class MapToSphere >
    struct GridObjectStreamTraits< SphereGrid< HostGrid, MapToSphere > >
    {
      typedef typename GridObjectStreamTraits< SphereGrid< HostGrid, MapToSphere > >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< SphereGrid< HostGrid, MapToSphere > >::OutStreamType OutStreamType;
    };
#endif // #if HAVE_DUNE_METAGRID

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDOBJECTSTREAMS_HH
