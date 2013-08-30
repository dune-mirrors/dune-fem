#ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BUTCHERTABLE_HH
#define DUNE_FEM_SOLVER_RUNGEKUTTA_BUTCHERTABLE_HH

//- dune-common includes 
#include <dune/common/dynmatrix.hh>
#include <dune/common/dynvector.hh>

namespace DuneODE 
{

  // SimpleButcherTable
  // ------------------

  template< class Field >
  class SimpleButcherTable
  {
    typedef SimpleButcherTable< Field > This;

  public:
    typedef Field FieldType;

    SimpleButcherTable ( int stages, int order, const FieldType *a, const FieldType *b, const FieldType *c )
    : stages_( stages ), order_( order ),
      a_( a ), b_( b ), c_( c )
    {}

    Dune::DynamicMatrix< FieldType > A () const { return makeMatrix( stages_, stages_, a_ ); }
    Dune::DynamicVector< FieldType > b () const { return makeVector( stages_, b_ ); }
    Dune::DynamicVector< FieldType > c () const { return makeVector( stages_, c_ ); }

    int order () const { return order_; }
    int stages () const { return stages_; }

  private:
    static Dune::DynamicMatrix< FieldType > makeMatrix ( int m, int n, const FieldType *data )
    {
      Dune::DynamicMatrix< FieldType > A( m, n );
      for( int i = 0; i < m; ++i )
        std::copy( data + i*n, data + (i+1)*n, A[ i ].begin() );
      return A;
    }

    static Dune::DynamicVector< FieldType > makeVector ( int n, const FieldType *data )
    {
      Dune::DynamicVector< FieldType > v( n );
      std::copy( data, data + n, v.begin() );
      return v;
    }

    int stages_, order_;
    const FieldType *a_, *b_, *c_;
  };



  // implicit butcher tables
  // -----------------------

  SimpleButcherTable< double > implicit34ButcherTable ();
  SimpleButcherTable< double > implicit3ButcherTable ();
  SimpleButcherTable< double > implicitEulerButcherTable ();
  SimpleButcherTable< double > gauss2ButcherTable ();



  // semiimplicit butcher tables
  // ---------------------------

  SimpleButcherTable< double > semiImplicitEulerButcherTable ( bool expl );
  SimpleButcherTable< double > semiImplicit23ButcherTable ( bool expl );
  SimpleButcherTable< double > semiImplicit33ButcherTable ( bool expl );
  SimpleButcherTable< double > semiImplicitSSP222ButcherTable ( bool expl );
  SimpleButcherTable< double > semiImplicitARK34ButcherTable ( bool expl );
  SimpleButcherTable< double > semiImplicitARK46ButcherTable ( bool expl );
  SimpleButcherTable< double > semiImplicitIERK45ButcherTable ( bool expl );

} // namespace DuneODE

#endif // #ifndef DUNE_FEM_SOLVER_RUNGEKUTTA_BUTCHERTABLE_HH
