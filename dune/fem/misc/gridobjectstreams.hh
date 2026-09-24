#ifndef DUNE_FEM_GRIDOBJECTSTREAMS_HH
#define DUNE_FEM_GRIDOBJECTSTREAMS_HH

#include <type_traits>
#include <dune/common/exceptions.hh>

#include <dune/fem/common/utility.hh>
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
      // generates a check for whether Grid has a type ObjectStreamType or not
      // the resulting class is hasMember_ObjectStreamType
      GENERATE_MEMBER_CHECK( ObjectStreamType );

      // for grids that export ObjectStreamType (e.g. ALUGrid etc)
      template <class G, bool >
      struct Selector
      {
        typedef typename G::ObjectStreamType type;
      };

      // for other grids that DO NOT export ObjectStreamType (e.g. YaspGrid etc)
      template <class G >
      struct Selector< G, false >
      {
        typedef DummyObjectStream type;
      };

      typedef typename Selector< Grid, hasMember_ObjectStreamType< Grid > :: value > :: type ObjectStreamType;

      typedef ObjectStreamType InStreamType;
      typedef ObjectStreamType OutStreamType;
    };

    // const version
    template< class Grid >
    struct GridObjectStreamTraits< const Grid >
    {
      typedef typename GridObjectStreamTraits< Grid >::InStreamType InStreamType;
      typedef typename GridObjectStreamTraits< Grid >::OutStreamType OutStreamType;
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_GRIDOBJECTSTREAMS_HH
