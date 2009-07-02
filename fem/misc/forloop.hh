#ifndef DUNE_FEM_FORLOOP_HH
#define DUNE_FEM_FORLOOP_HH

#include <dune/common/static_assert.hh>

namespace Dune
{

  template< template< int > class Operation, int first, int last >
  struct ForLoop
  {
    static void apply ()
    {
      Operation< first >::apply();
      ForLoop< Operation, first+1, last >::apply();
    }

    template< class T1 >
    static void apply ( T1 &p1 )
    {
      Operation< first >::apply( p1 );
      ForLoop< Operation, first+1, last >::apply( p1 );
    }
    
    template< class T1, class T2 >
    static void apply ( T1 &p1, T2 &p2 )
    {
      Operation< first >::apply( p1, p2 );
      ForLoop< Operation, first+1, last >::apply( p1, p2 );
    }

    template< class T1, class T2, class T3 >
    static void apply ( T1 &p1, T2 &p2, T3 &p3 )
    {
      Operation< first >::apply( p1, p2, p3 );
      ForLoop< Operation, first+1, last >::apply( p1, p2, p3 );
    }

    template< class T1, class T2, class T3, class T4 >
    static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4 )
    {
      Operation< first >::apply( p1, p2, p3, p4 );
      ForLoop< Operation, first+1, last >::apply( p1, p2, p3, p4 );
    }

    template< class T1, class T2, class T3, class T4, class T5 >
    static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4, T5 &p5 )
    {
      Operation< first >::apply( p1, p2, p3, p4, p5 );
      ForLoop< Operation, first+1, last >::apply( p1, p2, p3, p4, p5 );
    }

    template< class T1, class T2, class T3, class T4, class T5, class T6 >
    static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4, T5 &p5, T6 &p6 )
    {
      Operation< first >::apply( p1, p2, p3, p4, p5, p6 );
      ForLoop< Operation, first+1, last >::apply( p1, p2, p3, p4, p5, p6 );
    }

    template< class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
    static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4, T5 &p5, T6 &p6, T7 &p7 )
    {
      Operation< first >::apply( p1, p2, p3, p4, p5, p6, p7 );
      ForLoop< Operation, first+1, last >::apply( p1, p2, p3, p4, p5, p6, p7 );
    }

    template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
    static void apply ( T1 &p1, T2 &p2, T3 &p3, T4 &p4, T5 &p5, T6 &p6, T7 &p7, T8 &p8 )
    {
      Operation< first >::apply( p1, p2, p3, p4, p5, p6, p7, p8 );
      ForLoop< Operation, first+1, last >::apply( p1, p2, p3, p4, p5, p6, p7, p8 );
    }

  private:
    dune_static_assert( (first <= last), "ForLoop: first > last" );
  };


  template< template< int > class Operation, int last >
  struct ForLoop< Operation, last, last >
  {
    static void apply ()
    {
      Operation< last >::apply();
    }

    template< class T1 >
    static void apply ( T1 &p1 )
    {
      Operation< last >::apply( p1 );
    }
    
    template< class T1, class T2 >
    static void apply ( T1 &p1, T2 &p2 )
    {
      Operation< last >::apply( p1, p2 );
    }

    template< class T1, class T2, class T3 >
    static void apply ( T1 &p1, T2 &p2,  T3 &p3 )
    {
      Operation< last >::apply( p1, p2, p3 );
    }

    template< class T1, class T2, class T3, class T4 >
    static void apply ( T1 &p1, T2 &p2,  T3 &p3, T4 &p4 )
    {
      Operation< last >::apply( p1, p2, p3, p4 );
    }

    template< class T1, class T2, class T3, class T4, class T5 >
    static void apply ( T1 &p1, T2 &p2,  T3 &p3, T4 &p4, T5 &p5 )
    {
      Operation< last >::apply( p1, p2, p3, p4, p5 );
    }

    template< class T1, class T2, class T3, class T4, class T5, class T6 >
    static void apply ( T1 &p1, T2 &p2,  T3 &p3, T4 &p4, T5 &p5, T6 &p6 )
    {
      Operation< last >::apply( p1, p2, p3, p4, p5, p6 );
    }

    template< class T1, class T2, class T3, class T4, class T5, class T6, class T7 >
    static void apply ( T1 &p1, T2 &p2,  T3 &p3, T4 &p4, T5 &p5, T6 &p6, T7 &p7 )
    {
      Operation< last >::apply( p1, p2, p3, p4, p5, p6, p7 );
    }

    template< class T1, class T2, class T3, class T4, class T5, class T6, class T7, class T8 >
    static void apply ( T1 &p1, T2 &p2,  T3 &p3, T4 &p4, T5 &p5, T6 &p6, T7 &p7, T8 &p8 )
    {
      Operation< last >::apply( p1, p2, p3, p4, p5, p6, p7, p8 );
    }
  };

}

#endif
