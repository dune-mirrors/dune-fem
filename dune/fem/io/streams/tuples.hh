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
      InStreamFunctor functor;
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
      Dune::ForEachValue< Tuple > forEach( const_cast< Tuple >( tuple ) );
      OutStreamFunctor functor;
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



    // Dune::tuple to InStream
    // -----------------------

    template< class StreamTraits >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple<> &tuple )
    {
      return in;
    }

    template< class StreamTraits, class T1 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1 > &tuple )
    {
      typedef Dune::tuple< T1 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2 > &tuple )
    {
      typedef Dune::tuple< T1, T2 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2, T3 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2, T3, T4 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2, T3, T4, T5 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2, T3, T4, T5, T6 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2, T3, T4, T5, T6, T7 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6, T7 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9 >
    inline InStreamInterface< StreamTraits > &
      operator>> ( InStreamInterface< StreamTraits > &in, Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > Tuple;
      return TupleToInStream< InStreamInterface< StreamTraits > >::template apply< Tuple >( in, tuple );
    }



    // Dune::tuple to OutStream
    // ------------------------

    template< class StreamTraits >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple<> &tuple )
    {
      return out;
    }

    template< class StreamTraits, class T1 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1 > &tuple )
    {
      typedef Dune::tuple< T1 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2 > &tuple )
    {
      typedef Dune::tuple< T1, T2 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2, T3 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2, T3, T4 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2, T3, T4, T5 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2, T3, T4, T5, T6 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2, T3, T4, T5, T6, T7 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6, T7 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

    template< class StreamTraits, class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8, class T9 >
    inline OutStreamInterface< StreamTraits > &
      operator<< ( OutStreamInterface< StreamTraits > &out, const Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > &tuple )
    {
      typedef Dune::tuple< T1, T2, T3, T4, T5, T6, T7, T8, T9 > Tuple;
      return TupleToOutStream< OutStreamInterface< StreamTraits > >::template apply< Tuple >( out, tuple );
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_IO_STREAMS_TUPLES_HH
