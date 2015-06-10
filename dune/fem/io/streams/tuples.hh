#ifndef DUNE_FEM_IO_STREAMS_TUPLES_HH
#define DUNE_FEM_IO_STREAMS_TUPLES_HH

#include <dune/common/tuples.hh>
#include <dune/common/tupleutility.hh>

#include <dune/fem/io/streams/streams.hh>

namespace
{

  // TupleToInStream
  // ---------------

  template< class InStream >
  class TupleToInStream
  {
    struct InStreamFunctor
    {
      InStreamFunctor ( InStream &in ) : in_( in ) {}

      template< class T >
      void visit ( T &t ) const
      {
        in_ >> t;
      }

      private:
        InStream &in_;
    };

  public:
    template< class Tuple >
    static InStream &apply ( InStream &in, Tuple &tuple )
    {
      Dune::ForEachValue< Tuple > forEach( tuple );
      InStreamFunctor functor( in );
      forEach.apply( functor );
      return in;
    }
  };



  // TupleToOutStream
  // ----------------

  template< class OutStream >
  class TupleToOutStream
  {
    struct OutStreamFunctor
    {
      OutStreamFunctor ( OutStream &out ) : out_( out ) {}

      template< class T >
      void visit ( T &t ) const
      {
        out_ << t;
      }

      private:
        OutStream &out_;
    };

  public:
    template< class Tuple >
    static OutStream &apply ( OutStream &out, const Tuple &tuple )
    {
      Dune::ForEachValue< Tuple > forEach( const_cast< Tuple & >( tuple ) );
      OutStreamFunctor functor( out );
      forEach.apply( functor );
      return out;
    }
  };

} // namespace



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
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< std::tuple< Args... > >( in, tuple );
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
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< std::tuple < Args... > >( out, tuple );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_STREAMS_TUPLES_HH
