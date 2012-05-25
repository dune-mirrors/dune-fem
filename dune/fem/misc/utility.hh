#ifndef DUNE_FEM_UTILITY_HH
#define DUNE_FEM_UTILITY_HH

#include <dune/common/tupleutility.hh>

#include <dune/fem/misc/femtuples.hh>

namespace Dune
{

  /** @ addtogroup Common
   *
   * @{
   */

  /**
   * @file
   * @brief Contain utility classes which can be used with tuples.
   */


  /**
   * @brief Helper template to clone the type definition of a tuple with the
   * storage types replaced by a user-defined rule.
   *  
   * Suppose all storage types A_i in a tuple define a type A_i::B. You can
   * build up a pair consisting of the types defined by A_i::B in the following
   * way:
   \code
   template <class A>
   struct MyEvaluator {
     typedef typename A::B Type;
   };

   typedef ForEachType<MyEvaluator, ATuple>::Type BTuple;
   \endcode
   * Here, MyEvaluator is a helper struct that extracts the correct type from
   * the storage types of the tuple defined by the tuple ATuple.
   */

  template< template< class > class TypeEvaluator, class TupleType >
  struct ForEachTupleType 
  /** \cond */
  : ForEachTupleType< TypeEvaluator, typename TupleType :: FirstPair >
  /** \endcond */
  {};

  // Specialisation for standard tuple element
  template <template <class> class TypeEvaluator, class Head, class Tail>
  struct ForEachTupleType<TypeEvaluator, Pair<Head, Tail> > {
    //! Defines type corresponding to the subtuple defined by Pair<Head, Tail>
    typedef Pair<typename TypeEvaluator<Head>::Type, 
                 typename ForEachTupleType<TypeEvaluator, Tail>::Type> Type;
  };
  // Specialisation for last element
  template <template <class> class TypeEvaluator>
  struct ForEachTupleType<TypeEvaluator, Nil> {
    typedef Nil Type;
  };

//   template <template <class> class TypeEvaluator, class Head>
//   struct ForEachType<TypeEvaluator, Pair<Head, Nil> > {
//     //! For the storage element, Head is replaced by the expression provided
//     //! by the TypeEvaluator helper template.
//     typedef Pair<typename TypeEvaluator<Head>::Type, Nil> Type;
//   };

  /**
   * @brief Helper template which implements iteration over all storage 
   * elements in a tuple.
   *
   * Compile-time constructs that allows to process all elements in a tuple.
   * The exact operation performed on an element is defined by a function
   * object, which needs to implement a visit method which is applicable to
   * all storage elements of a tuple.
   *
   * The following example implements a function object which counts the
   * elements in a tuple
   \code
   template <class T>
   struct Counter {
   Counter() : result_(0) {}
   
   template <class T>
   void visit(T& elem) { ++result_; }

   int result_;
   };
   \endcode
   * The number of elements in the tuple are stored in the member variable 
   * result_. The Counter can be used as follows, assuming a tuple t of type
   * MyTuple is given:
   \code
   Counter c;
   ForEachValue<MyTuple> forEach(t);

   forEach.apply(c);
   std::cout << "Number of elements is: " << c.result_ << std::endl;
   \endcode
   */

  template <class TupleType>
  class ForEachTupleValue
  {
  public:
    //! \brief Constructor
    //! \param tuple The tuple which we want to process.
    ForEachTupleValue(TupleType& tuple) : tuple_(tuple) {}
    
    //! \brief Applies a function object to each storage element of the tuple.
    //! \param f Function object.
    template <class Functor>
    void apply(Functor& f) {
      apply(f, tuple_);
    }
    
  private:
    //! Specialisation for the last element
    template <class Functor, class Head>
    void apply(Functor& f, Pair<Head, Nil>& last) {
      f.visit(last.first());
    }
    
    //! Specialisation for a standard tuple element
    template <class Functor, class Head, class Tail>
    void apply(Functor& f, Pair<Head, Tail>& pair) {
      f.visit(pair.first());
      apply(f, pair.second());
    }
  private:
    TupleType& tuple_;
  };

  //- Definition ForEachValuePair class
  // Assertion: both tuples have the same length and the contained types are
  // compatible in the sense of the applied function object
  /**
   * @brief Extension of ForEachValue to two tuples...
   *
   * This class provides the framework to process two tuples at once. It works
   * the same as ForEachValue, just that the corresponding function object
   * takes one argument from the first tuple and one argument from the second.
   *
   * \note You have to ensure that the two tuples you provide are compatible
   * in the sense that they have the same length and that the objects passed
   * to the function objects are related in meaningful way. The best way to
   * enforce it is to build the second tuple from the existing first tuple
   * using ForEachType.
   */
  template <class TupleType1, class TupleType2>
  class ForEachTupleValuePair
  {
  public:
    //! Constructor
    //! \param t1 First tuple.
    //! \param t2 Second tuple.
    ForEachTupleValuePair(TupleType1& t1, TupleType2& t2) :
      tuple1_(t1),
      tuple2_(t2)
    {}

    //! Applies the function object f to the pair of tuples.
    //! \param f The function object to apply on the pair of tuples.
    template <class Functor>
    void apply(Functor& f) {
      apply(f, tuple1_, tuple2_);
    }

  private:
    //! Specialisation for the last element.
    template <class Functor, class Head1, class Head2>
    void apply(Functor& f, Pair<Head1, Nil>& last1, Pair<Head2, Nil>& last2) {
      f.visit(last1.first(), last2.first());
    }

    //! Specialisation for a standard element.
    template <class Functor, class Head1, class Tail1, class Head2,class Tail2>
    void apply(Functor& f, Pair<Head1, Tail1>& p1, Pair<Head2, Tail2>& p2) {
      f.visit(p1.first(), p2.first());
      apply(f, p1.second(), p2.second());
    }

  private:
    TupleType1& tuple1_;
    TupleType2& tuple2_;
  };

  template <class T1,class T2>
  struct CombinePairs {
  };


  template <class TupType>
  struct TupleToPair {
    typedef typename TupleToPair<typename TupType::Type2>::ReturnType NextType;
    typedef Pair<typename TupType::Type1,NextType> ReturnType;
    static ReturnType convert(TupType tup) {
      NextType next = TupleToPair<typename TupType::Type2>::convert(tup.second());
      return ReturnType(tup.first(),next);
    }
  };
  template <class T>
  struct TupleToPair<Pair<T,Nil> > {
    typedef Pair<T,Nil> ReturnType;
    static ReturnType convert(Pair<T,Nil> tup) {
      return ReturnType(tup.first(),nullType());
    }
  };


}

#endif
