#ifndef DUNE_FEM_METAPROGRAMMING_HH
#define DUNE_FEM_METAPROGRAMMING_HH

namespace Dune
{

  namespace Fem
  {

    // meta programming for boolean values

    template< bool b >
    struct MetaBool
    {
      // for compatibility other dune modules, v is preferred
      static const bool v = b;
      static const bool value = v;
    };

    struct MetaTrue
    : public MetaBool< true >
    {};

    struct MetaFalse
    : public MetaBool< false >
    {};



    template< class A >
    struct MetaNot
    : public MetaBool< !(A :: v) >
    {};

    template< class A, class B >
    struct MetaAnd
    : public MetaBool< (A :: v && B :: v) >
    {};

    template< class A, class B >
    struct MetaOr
    : public MetaBool< (A :: v || B :: v) >
    {};



    template< bool condition, class Then, class Else >
    struct MetaIf;

    template< class Then, class Else >
    struct MetaIf< true, Then, Else >
    : public Then
    {};

    template< class Then, class Else >
    struct MetaIf< false, Then, Else >
    : public Else
    {};



    template< bool condition,
              template< unsigned int > class Then, unsigned int i,
              template< unsigned int > class Else, unsigned int j >
    struct MetaIfTemplate;

    template< template< unsigned int > class Then, unsigned int i,
              template< unsigned int > class Else, unsigned int j >
    struct MetaIfTemplate< true, Then, i, Else, j >
    : public Then< i >
    {};

    template< template< unsigned int > class Then, unsigned int i,
              template< unsigned int > class Else, unsigned int j >
    struct MetaIfTemplate< false, Then, i, Else, j >
    : public Else< j >
    {};



    template< bool condition >
    struct MetaAssert;

    template<>
    struct MetaAssert< true >
    {};



    // meta programming for integer values

    template< int i >
    struct MetaInt
    {
      static const int value = i;
      static const int v = i;
    };




    template< class A, class B >
    struct MetaPlus
    : public MetaInt< (A :: value + B :: value) >
    {};

    template< class A, class B >
    struct MetaMinus
    : public MetaInt< (A :: value - B :: value) >
    {};



    template< class A >
    class MetaInc
    : public MetaInt< (A :: value + 1) >
    {};

    template< class A >
    class MetaDec
    : public MetaInt< (A :: value - 1) >
    {};



    template< class A, class B >
    struct MetaMultiply
    : public MetaInt< (A :: value * B :: value) >
    {};

    template< class A, class B >
    struct MetaDivide
    : public MetaInt< (A :: value / B :: value) >
    {};



    template< class A, class B >
    struct MetaMin
    : public MetaInt< (A :: value <= B :: value ? A :: value : B :: value) >
    {};

    template< class A, class B >
    struct MetaMax
    : public MetaInt< (A :: value >= B :: value ? A :: value : B :: value) >
    {};



    template< class A, class B >
    struct MetaLess
    : public MetaBool< (A :: value < B :: value) >
    {};

    template< class A, class B >
    struct MetaLessEqual
    : public MetaBool< (A :: value <= B :: value) >
    {};

    template< class A, class B >
    struct MetaGreaterEqual
    : public MetaBool< (A :: value >= B :: value) >
    {};

    template< class A, class B >
    struct MetaGreater
    : public MetaBool< (A :: value > B :: value) >
    {};



    // Meta programming for functions

    template< class A, class B >
    struct MetaSequence
    {
      typedef typename A :: ArgumentType ArgumentType;

      inline static void apply ( ArgumentType argument )
      {
        A :: apply( argument );
        B :: apply( argument );
      }
    };



    // Loop: apply an operation i times (which means there are i+1 arguments)

    template< template< class, class > class Operation,
              template< unsigned int > class f,
              unsigned int i >
    struct Loop
    : public Operation< Loop< Operation, f, i-1 >, f< i > >
    {
    };

    template< template< class, class > class Operation,
              template< unsigned int > class f >
    struct Loop< Operation, f, 0 >
    : public f< 0 >
    {
    };



    // MetaProtect: protect template instanciation

    template< class Void,
              bool condition,
              template< unsigned int > class Value,
              unsigned int i >
    struct MetaProtect;

    template< class Void, template< unsigned int > class Value, unsigned int i >
    struct MetaProtect< Void, true, Value, i >
    : public Value< i >
    {
    };

    template< class Void, template< unsigned int > class Value, unsigned int i >
    struct MetaProtect< Void, false, Value, i >
    : public Void
    {
    };



    template< bool condition, template< unsigned int > class Value, unsigned int i >
    struct MetaProtectInt
    : public MetaProtect< MetaInt< 0 >, condition, Value, i >
    {
    };




    template< template< unsigned int > class A, unsigned int i,
              class B, unsigned int n >
    struct Protect
    : public A< i >
    {
    };

    template< template< unsigned int > class A, unsigned int i, class B >
    struct Protect< A, i, B, i >
    : public B
    {
    };

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_METAPROGRAMMING_HH
