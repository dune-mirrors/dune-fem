#ifndef DUNE_FEM_METAPROGRAMMING_HH
#define DUNE_FEM_METAPROGRAMMING_HH

namespace Dune
{
  
   // meta programming for boolean values
  
  template< bool b >
  struct MetaBool
  {
    enum { value = b };
  };



  template< class A >
  struct MetaNot
  : public MetaBool< !A :: value >
  {
  };



  template< class A, class B >
  struct MetaAnd
  : public MetaBool< A :: value && B :: value >
  {
  };
  
  template< class A, class B >
  struct MetaOr
  : public MetaBool< A :: value || B :: value >
  {
  };


 
  template< bool Condition, class Then, class Else >
  struct If;

  template< class Then, class Else >
  struct If< true, Then, Else >
  : public Then
  {
  };

  template< class Then, class Else >
  struct If< false, Then, Else >
  : public Else
  {
  };
 
  template< class Condition, class Then, class Else >
  struct MetaIf
  : public If< Condition :: value, Then, Else >
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
 
}

#endif
