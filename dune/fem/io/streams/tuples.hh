#ifndef DUNE_FEM_IO_STREAMS_TUPLES_HH
#define DUNE_FEM_IO_STREAMS_TUPLES_HH

#include <tuple>
#include <utility>

#include <dune/common/hybridutilities.hh>
#include <dune/fem/io/streams/streams.hh>


namespace Dune
{

  namespace Fem
  {

    // External forward declarations
    // -----------------------------

    template <class StreamTraits>
    class OutStreamInterface;
    template <class StreamTraits>
    class InStreamInterface;



    // std::tuple to InStream
    // ----------------------

    template< class StreamTraits >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, std::tuple<> &tuple )
    {
      return in;
    }

    template< class StreamTraits, class... Args >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, std::tuple< Args... > &tuple )
    {

      Hybrid::forEach ( std::make_index_sequence< sizeof...( Args ) >{}, [ & ]( auto i ){ in >> std::get< i >( tuple ); } );
      return in;
    }



    // std::tuple to OutStream
    // -----------------------

    template< class StreamTraits >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const std::tuple<> &tuple )
    {
      return out;
    }

    template< class StreamTraits, class... Args >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const std::tuple< Args... > &tuple )
    {
      Hybrid::forEach ( std::make_index_sequence< sizeof...( Args ) >{}, [ & ]( auto i ){ out << std::get< i >( tuple ); } );
      return out;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_STREAMS_TUPLES_HH
