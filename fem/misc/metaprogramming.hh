#ifndef DUNE_FEM_METAPROGRAMMING_HH
#define DUNE_FEM_METAPROGRAMMING_HH

namespace Dune
{
  
  // meta programming for boolean values
  
  struct MetaTrue
  {
    typedef MetaTrue Value;
  };

  struct MetaFalse
  {
    typedef MetaFalse Value;
  };



  template< class A >
  struct MetaNot
  : public MetaNot< typename A :: Value >
  {
  };
 
  template<>
  struct MetaNot< MetaTrue >
  : public MetaFalse
  {
  };

  template<>
  struct MetaNot< MetaFalse >
  : public MetaTrue
  {
  };



  template< class A, class B >
  struct MetaAnd
  : public MetaAnd< typename A :: Value, typename B :: Value >
  {
  };

  template<>
  struct MetaAnd< MetaTrue, MetaTrue >
  : public MetaTrue
  {
  };

  template<>
  struct MetaAnd< MetaTrue, MetaFalse >
  : public MetaFalse
  {
  };

  template<>
  struct MetaAnd< MetaFalse, MetaTrue >
  : public MetaFalse
  {
  };

  template<>
  struct MetaAnd< MetaFalse, MetaFalse >
  : public MetaFalse
  {
  };



  template< class A, class B >
  struct MetaOr
  : public MetaOr< typename A :: Value, typename B :: Value >
  {
  };

  template<>
  struct MetaOr< MetaTrue, MetaTrue >
  : public MetaTrue
  {
  };

  template<>
  struct MetaOr< MetaTrue, MetaFalse >
  : public MetaTrue
  {
  };

  template<>
  struct MetaOr< MetaFalse, MetaTrue >
  : public MetaTrue
  {
  };

  template<>
  struct MetaOr< MetaFalse, MetaFalse >
  : public MetaFalse
  {
  };



  template< bool >
  struct MetaBool;

  template<>
  struct MetaBool< true >
  : public MetaTrue
  {
  };

  template<>
  struct MetaBool< false >
  : public MetaFalse
  {
  };



  template< class Condition, class Then, class Else >
  struct MetaIf
  : public MetaIf< typename Condition :: Value, Then, Else >
  {
  };

  template< class Then, class Else >
  struct MetaIf< MetaTrue, Then, Else >
  : public Then
  {
  };

  template< class Then, class Else >
  struct MetaIf< MetaFalse, Then, Else >
  : public Else
  {
  };
 


  // meta programming for integer values
  
  template< int i >
  struct MetaInt
  {
    enum { value = i };
  };



 
  template< class A, class B >
  struct MetaPlus
  : public MetaInt< A :: value + B :: value >
  {
  };

  template< class A, class B >
  struct MetaMinus
  : public MetaInt< A :: value - B :: value >
  {
  };



  template< class A >
  class MetaInc
  : public MetaInt< A :: value + 1 >
  {
  };

  template< class A >
  class MetaDec
  : public MetaInt< A :: value - 1 >
  {
  };


  
  template< class A, class B >
  struct MetaMultiply
  : public MetaInt< A :: value * B :: value >
  {
  };

  template< class A, class B >
  struct MetaDivide
  : public MetaInt< A :: value / B :: value >
  {
  };

  
  
  template< class A, class B >
  struct MetaMin
  : public MetaInt< (A :: value <= B :: value ? A :: value : B :: value) >
  {
  };

  template< class A, class B >
  struct MetaMax
  : public MetaInt< (A :: value >= B :: value ? A :: value : B :: value) >
  {
  };



  template< class A, class B >
  struct MetaLess
  : public MetaBool< (A :: value < B :: value) >
  {
  };

  template< class A, class B >
  struct MetaLessEqual
  : public MetaBool< (A :: value <= B :: value) >
  {
  };
  
  template< class A, class B >
  struct MetaGreaterEqual
  : public MetaBool< (A :: value >= B :: value) >
  {
  };
  
  template< class A, class B >
  struct MetaGreater
  : public MetaBool< (A :: value > B :: value) >
  {
  };



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



  // Protect: protect template instanciation
  
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

}

#endif
