#ifndef DUNE_FEM_MISC_UMFPACK_HH
#define DUNE_FEM_MISC_UMFPACK_HH

#include <complex>

#if HAVE_DUNE_ISTL
#include <dune/istl/umfpack.hh>
#endif

#if HAVE_SUITESPARSE_UMFPACK
#include <umfpack.h>

#include <dune/fem/misc/double.hh>

namespace Dune
{

  // UMFPackMethodChooser
  // --------------------

#if not HAVE_DUNE_ISTL
  template< class T >
  struct UMFPackMethodChooser {};


  template< >
  struct UMFPackMethodChooser< double >
  {
    template< class... A >
    static void defaults ( A ... args )
    {
      umfpack_dl_defaults( args ... );
    }
    template< class ... A >
    static void free_numeric ( A ... args )
    {
      umfpack_dl_free_numeric( args ... );
    }
    template< class ... A >
    static void free_symbolic ( A ... args )
    {
      umfpack_dl_free_symbolic( args ... );
    }
    template< class ... A >
    static int load_numeric ( A ... args )
    {
      return umfpack_dl_load_numeric( args ... );
    }
    template< class ... A >
    static void numeric ( A ... args )
    {
      umfpack_dl_numeric( args ... );
    }
    template< class ... A >
    static void report_info ( A ... args )
    {
      umfpack_dl_report_info( args ... );
    }
    template< class ... A >
    static void report_status ( A ... args )
    {
      umfpack_dl_report_status( args ... );
    }
    template< class ... A >
    static int save_numeric ( A ... args )
    {
      return umfpack_di_save_numeric( args ... );
    }
    template< class ... A >
    static void solve ( A ... args )
    {
      umfpack_dl_solve( args ... );
    }
    template< class ... A >
    static void symbolic ( A ... args )
    {
      umfpack_dl_symbolic( args ... );
    }
  };

  template< >
  struct UMFPackMethodChooser< std::complex< double > >
  {
    template< class ... A >
    static void defaults ( A ... args )
    {
      umfpack_zi_defaults( args ... );
    }
    template< class ... A >
    static void free_numeric ( A ... args )
    {
      umfpack_zi_free_numeric( args ... );
    }
    template< class ... A >
    static void free_symbolic ( A ... args )
    {
      umfpack_zi_free_symbolic( args ... );
    }
    template< class ... A >
    static int load_numeric ( A ... args )
    {
      return umfpack_zi_load_numeric( args ... );
    }
    template< class ... A >
    static void numeric ( const int *cs, const int *ri, const double *val, A ... args )
    {
      umfpack_zi_numeric( cs, ri, val, nullptr, args ... );
    }
    template< class ... A >
    static void report_info ( A ... args )
    {
      umfpack_zi_report_info( args ... );
    }
    template< class ... A >
    static void report_status ( A ... args )
    {
      umfpack_zi_report_status( args ... );
    }
    template< class ... A >
    static int save_numeric ( A ... args )
    {
      return umfpack_zi_save_numeric( args ... );
    }
    template< class ... A >
    static void solve ( int m, const int *cs, const int *ri, std::complex< double > *val, double *x, const double *b, A ... args )
    {
      const double *cval = reinterpret_cast< const double * >(val);
      umfpack_zi_solve( m, cs, ri, cval, nullptr, x, nullptr, b, nullptr, args ... );
    }
    template< class ... A >
    static void symbolic ( int m, int n, const int *cs, const int *ri, const double *val, A ... args )
    {
      umfpack_zi_symbolic( m, n, cs, ri, val, nullptr, args ... );
    }
  };
#endif // #if not HAVE_DUNE_ISTL

  template< >
  struct UMFPackMethodChooser< Fem::Double >
  {
    template< class... A >
    static void defaults ( A ... args )
    {
      umfpack_di_defaults( args ... );
    }
    template< class ... A >
    static void free_numeric ( A ... args )
    {
      umfpack_di_free_numeric( args ... );
    }
    template< class ... A >
    static void free_symbolic ( A ... args )
    {
      umfpack_di_free_symbolic( args ... );
    }
    template< class ... A >
    static int load_numeric ( A ... args )
    {
      return umfpack_di_load_numeric( args ... );
    }
    template< class ... A >
    static void numeric ( A ... args )
    {
      umfpack_di_numeric( args ... );
    }
    template< class ... A >
    static void report_info ( A ... args )
    {
      umfpack_di_report_info( args ... );
    }
    template< class ... A >
    static void report_status ( A ... args )
    {
      umfpack_di_report_status( args ... );
    }
    template< class ... A >
    static int save_numeric ( A ... args )
    {
      return umfpack_di_save_numeric( args ... );
    }
    template< class ... A >
    static void symbolic ( A ... args )
    {
      umfpack_di_symbolic( args ... );
    }
    static void numeric ( const int *Ap, const int *Ai, const double *Ax, void* Symbolic, void **Numeric,
        const double Control[ UMFPACK_CONTROL], double Info[ UMFPACK_INFO] )
    {
      umfpack_di_numeric( Ap, Ai, Ax, Symbolic, Numeric, Control, Info );
      Fem::FlOpCounter< Fem::Double >::instance() += Info[ UMFPACK_FLOPS ];
    }

    static void solve ( int m, const int *cs, const int *ri, const Fem::Double *val, Fem::Double *x, const Fem::Double *b,
        void *Numeric, const double Control[ UMFPACK_CONTROL ], double Info[ UMFPACK_INFO ] )
    {
      umfpack_di_solve( m, cs, ri, reinterpret_cast< const double * >(val),
          reinterpret_cast< double * >( x ), reinterpret_cast< const double * >( b ), Numeric, Control, Info );
      Fem::FlOpCounter< Fem::Double >::instance() += Info[ UMFPACK_SOLVE_FLOPS ];
    }
  };

} // namespace Dune

#endif // #if HAVE_SUITESPARSE_UMFPACK

#endif // #ifndef DUNE_FEM_MISC_UMFPACK_HH
