#ifndef DUNE_CALLER_HH
#define DUNE_CALLER_HH

#include <dune/common/tuples.hh>
#include <dune/common/utility.hh>
#include <dune/common/typetraits.hh>

namespace Dune {

  // Template magic abound here...
  /**
   * @brief Extracts elements from a tuple based on a Selector definition
   *
   * Use this class to extract the elements from a tuple. The selector 
   * specifies the elements to be extracted by indicating their index. The
   * extracted elements are again stored in a tuple.
   */
  template <class ArgTupleImp, class SelectorImp>
  struct Filter 
  {
    typedef SelectorImp SelectorType;
    //! The index of the extracted element.
    enum { index = ElementType<0, SelectorType>::Type::value };
    //! The type of the extracted element.
    typedef typename ElementType<index, ArgTupleImp>::Type AppendType;
    //! The type of the next stage of the filter.
    typedef Filter<ArgTupleImp, typename SelectorType::Type2> NextFilterType;
    //! The filtered tuple resulting from this filter stage.
    typedef Pair<AppendType, typename NextFilterType::ResultType> ResultType;
    
    //! Extracts the elements specified by the selector template argument.
    //! \param arg Argument tuple.
    static inline ResultType apply(ArgTupleImp& arg) {
      return ResultType(Element<index>::get(arg), NextFilterType::apply(arg));
    }
  };
  
  /**
   * @brief Specialisation for closure.
   *
   * The end of the filtering is marked by the last element of the selector.
   * \warning If the largest index in the selector is greater than the length
   * of the tuple, an ugly compile-time error will occur... (I'll switch that
   * off as soon as I can)
   */ 
  template <class ArgTupleImp>
  struct Filter<ArgTupleImp, Nil>
  {
    //! The closure of a tuple is always Nil
    typedef Nil ResultType;
    
    //! Provides the closure of the filtered tuple
    static inline ResultType apply(ArgTupleImp& arg) {
      return Nil();
    }
  };

  //- Helper class for helper class below ;-)
  // ! Need to strip pointer first, then add it again

  /**
   * @brief Extracts the type of the LocalFunction pointers from a tuple of 
   * DiscreteFunction pointers.
   */
  template <class DFType>
  struct DFTypeEvaluator {
    typedef typename TypeTraits<DFType>::PointeeType::LocalFunctionType Type;
  };

  /**
   * @brief Extracts the type of the Range vector pointers from a tuple of 
   * LocalFunction pointers.
   */
  template <class LFType>
  struct RangeTypeEvaluator {
    typedef typename LFType::RangeType Type;
  };


  // * Idea: can I do the same with ForEachValue?
  // * Give tuple in constructor and get the rest with apply()
  /*
  struct LocalFunctionCreator2 {
    template <class Head, class Tail>
    struct Traits {
      typedef typename ForEachType<DFTypeEvaluator, Tail>::Type TailType;
      typedef typename ForEachType<
        DFTypeEvaluator, Pair<Head, Tail> >::Type ResultType;
      typedef typename TypeTraits<
        Head>::PointeeType::LocalFunctionType LocalFunctionType;
    };

    template <class Head>
    static inline typename Traits<Head, Nil>::ResultType
    apply(Pair<Head, Nil>& p) {
      typedef typename Traits<Head, Nil>::ResultType Result;
      typedef typename Traits<Head, Nil>::LocalFunctionType LocalFunction;
      return Result(new LocalFunction(p.first()->newLocalFunction()),
                    Nil());
    }

    template <class Head, class Tail>
    static inline typename Traits<Head, Tail>::ResultType
    apply(Pair<Head, Tail>& p) {
      typedef typename Traits<Head, Tail>::ResultType Result;
      typedef typename Traits<Head, Tail>::LocalFunctionType LocalFunction;
      return Result(new LocalFunction(p.first()->newLocalFunction()), 
                    LocalFunctionCreator2::apply(p.second()));
    }
  };
  */
  /**
   * @brief Creates a tuple of local function pointers from a tuple of 
   * discrete function pointers by calling newLocalFunction on the discrete
   * function (plus a copy constructor...)
   * \note Don't forget to delete the local functions in the end!
   */
  template <class DFTupleType>
  class LocalFunctionCreator {};

  /**
   * @brief Specialisation for standard tuple element
   */
  template <class Head, class Tail>
  class LocalFunctionCreator<Pair<Head, Tail> > {
  private: 
    template <class T>
    friend class LocalFunctionCreator;
    typedef typename ForEachType<DFTypeEvaluator, Tail>::Type TailType;
  
  public:
    typedef typename ForEachType<
      DFTypeEvaluator, Pair<Head, Tail> >::Type ResultType;
    typedef typename TypeTraits<
      Head>::PointeeType::LocalFunctionType LocalFunctionType;

  public:
    static inline ResultType apply(Pair<Head, Tail>& pairs) {
      return 
        ResultType(LocalFunctionType(*pairs.first()),
                   LocalFunctionCreator<Tail>::apply(pairs.second()));
    }
  };

  /**
   * @brief Specialisation for last element
   */
  template <class Head>
  class LocalFunctionCreator<Pair<Head, Nil> > {
  private:
    template <class T>
    friend class LocalFunctionCreator;
  
  public:
    typedef typename TypeTraits<Head>::PointeeType::LocalFunctionType LocalFunctionType;
    typedef Pair<LocalFunctionType, Nil> ResultType;
 
  public:
    static inline ResultType apply(Pair<Head, Nil>& pairs) {
      return 
        ResultType(LocalFunctionType(*pairs.first()),
                   Nil());
    }
  };

  /**
   * @brief Generates a tuple with range vectors
   */
  template <class LFTupleType>
  class RangeVectorCreator {};
  
  /**
   * @brief Specialisation for standard element
   */
  template <class Head, class Tail>
  class RangeVectorCreator<Pair<Head, Tail> > {
  private:
    template <class T>
    friend class RangeVectorCreator;
    typedef typename ForEachType<RangeTypeEvaluator, Tail>::Type TailType;
  
  public:
    typedef typename ForEachType<
      RangeTypeEvaluator, Pair<Head, Tail> >::Type ResultType;
    typedef typename Head::RangeType RangeType;

  public:
    static inline ResultType apply() {
      return ResultType(RangeType(0.0), 
                        RangeVectorCreator<Tail>::apply());
    }

    static inline ResultType apply(Pair<Head, Tail>& pairs) {
      return ResultType(RangeType(0.0), 
                        RangeVectorCreator<Tail>::apply(pairs.second()));
    }
  };
  
  /**
   * @brief Specialisation for last element
   */
  template <class Head>
  class RangeVectorCreator<Pair<Head, Nil> > {
    template <class T>
    friend class RangeVectorCreator;

  public:
    typedef typename Head::RangeType RangeType;
    typedef Pair<RangeType, Nil> ResultType;

  public:
    static inline ResultType apply() {
      return ResultType(RangeType(0.0), Nil());
    }

    static inline ResultType apply(Pair<Head, Nil>& pairs) {
      return ResultType(RangeType(0.0), Nil());
    }
  };

  /**
   * @brief Functor class which sets a tuple of local functions to a specific
   * entity.
   *
   * This Functor as an argument to ForEachValuePair with a tuple of discrete
   * function pointers and the corresponding local function pointers.
   */
  template <class EntityImp>
  class LocalFunctionSetter {
  public:
    //! Constructor
    //! \param en The entity the local functions are set to.
    LocalFunctionSetter(EntityImp& en) :
      en_(en)
    {}

    //! Applies the setting on every DiscreteFunction/LocalFunction pair.
    template <class DFPointer, class LFType>
    void visit(DFPointer df, LFType& lf) {
      lf = df->localFunction(en_);
    }

  private:
    EntityImp& en_;
  };

  /**
   * @brief Calls lf.evaluateLocal on a tuple of local functions with a 
   * corresponding tuple of results.
   *
   * Use this Functor in conjunction with a ForEachValuePair built from a tuple
   * of local functions and a tuple with the corresponding range vectors.
   */
  template <class EntityImp, class DomainImp>
  class LocalFunctionEvaluateLocal {
  public:
    //! Constructor
    //! \param en Entity on which the local function is evaluated.
    //! \param x The local coordinate on en.
    LocalFunctionEvaluateLocal(EntityImp& en, DomainImp& x) :
      en_(en),
      x_(x)
    {}

    //! Triggers the evaluation of a local function
    template <class LFType, class RangeType>
    void visit(LFType& lf, RangeType& res) {
      lf.evaluateLocal(en_, x_, res);
    }

  private:
    EntityImp& en_;
    DomainImp& x_;
  };

  /**
   * @brief  Calls lf.evaluate on a tuple of local functions with a 
   * corresponding tuple of results.
   *
   * Use this Functor in conjunction with a ForEachValuePair built from a tuple
   * of local functions and a tuple with the corresponding range vectors.
   *
   */
  template <class EntityImp, class QuadratureImp>
  class LocalFunctionEvaluateQuad {
  public:
    //! Constructor
    //! \param en The entity the local functions are evaluated on.
    //! \param quad The quadrature in question.
    //! \param quadPoint The index of the quadrature point of quadrature quad
    LocalFunctionEvaluateQuad(EntityImp& en, 
                              QuadratureImp& quad,
                              int quadPoint) :
      en_(en),
      quad_(quad),
      quadPoint_(quadPoint)
    {}

    //! Evaluation of a local function
    template <class LFType, class RangeType>
    void visit(LFType& lf, RangeType& res) {
      lf.evaluate(en_, quad_, quadPoint_, res);
    }

  private:
    EntityImp& en_;
    QuadratureImp& quad_;
    int quadPoint_;
  };

  /**
   * @brief Helper class which actually calls the functions of the problem
   *
   * This class enables us to call methods by varying argument lengths.
   */
  /*
  template <class RangePairImp>
  class Caller {};
  
  template <class T1>
  class Caller<Pair<T1, Nil> > {
  public:
    typedef Pair<T1, Nil> RangePairType;
    
    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               const DomainType& x,
                               RangePairType& p, 
                               ResultType& res) {
      problem.analyticalFlux(en, x, 
                             Element<0>::get(p), 
                             res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               QuadratureType& quad,
                               int quadPoint,
                               RangePairType& p,
                               ResultType& res) {
      problem.analyticalFlux(en, quad, quadPoint,
                             Element<0>::get(p), 
                             res);
    }

    template <
      class ProblemType, 
      class IntersectionIterator, 
      class DomainType, 
      class ResultType
      >
    static void numericalFluxFlux(ProblemType& problem,
                                  IntersectionIterator& nit,
                                  const DomainType& x,
                                  RangePairType& pEn,
                                  RangePairType& pNeigh,
                                  ResultType& resEn, 
                                  ResultType& resNeigh) {
      problem.numericalFlux(nit, x, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType,
      class IntersectionIterator,
      class QuadratureType,
      class ResultType
      >
    static void numericalFlux(ProblemType& problem,
                              IntersectionIterator& nit,
                              QuadratureType& quad,
                              int quadPoint,
                              RangePairType& pEn,
                              RangePairType& pNeigh,
                              ResultType& resEn,
                              ResultType& resNeigh) {
      problem.numericalFlux(nit, quad, quadPoint, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       const DomainType& x,
                       RangePairType& p,
                       ResultType& res) {
      problem.source(en, x, 
                     Element<0>::get(p), 
                     res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       QuadratureType& quad,
                       int quadPoint,
                       RangePairType& p,
                       ResultType res) {
      problem.source(en, quad, quadPoint, 
                     Element<0>::get(p), 
                     res);
    }
  };
  
  // Two arguments
  template <class T1, class T2>
  class Caller<Pair<T1, Pair<T2, Nil> > > {
  public:
    typedef Pair<T1, Pair<T2, Nil> > RangePairType;
    
    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               const DomainType& x,
                               RangePairType& p, 
                               ResultType& res) {
      problem.analyticalFlux(en, x, 
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               QuadratureType& quad,
                               int quadPoint,
                               RangePairType& p,
                               ResultType& res) {
      problem.analyticalFlux(en, quad, quadPoint,
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             res);
    }

    template <
      class ProblemType, 
      class IntersectionIterator, 
      class DomainType, 
      class ResultType
      >
    static void numericalFluxFlux(ProblemType& problem,
                                  IntersectionIterator& nit,
                                  const DomainType& x,
                                  RangePairType& pEn,
                                  RangePairType& pNeigh,
                                  ResultType& resEn, 
                                  ResultType& resNeigh) {
      problem.numericalFlux(nit, x, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType,
      class IntersectionIterator,
      class QuadratureType,
      class ResultType
      >
    static void numericalFlux(ProblemType& problem,
                              IntersectionIterator& nit,
                              QuadratureType& quad,
                              int quadPoint,
                              RangePairType& pEn,
                              RangePairType& pNeigh,
                              ResultType& resEn,
                              ResultType& resNeigh) {
      problem.numericalFlux(nit, quad, quadPoint, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       const DomainType& x,
                       RangePairType& p,
                       ResultType& res) {
      problem.source(en, x, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       QuadratureType& quad,
                       int quadPoint,
                       RangePairType& p,
                       ResultType res) {
      problem.source(en, quad, quadPoint, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     res);
    }
  };

  // Three arguments
  template <class T1, class T2, class T3>
  class Caller<Pair<T1, Pair<T2, Pair<T3, Nil> > > > {
  public:
    typedef Pair<T1, Pair<T2, Pair<T3, Nil> > > RangePairType;
    
    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               const DomainType& x,
                               RangePairType& p, 
                               ResultType& res) {
      problem.analyticalFlux(en, x, 
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             Element<2>::get(p), 
                             res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               QuadratureType& quad,
                               int quadPoint,
                               RangePairType& p,
                               ResultType& res) {
      problem.analyticalFlux(en, quad, quadPoint,
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             Element<2>::get(p), 
                             res);
    }

    template <
      class ProblemType, 
      class IntersectionIterator, 
      class DomainType, 
      class ResultType
      >
    static void numericalFluxFlux(ProblemType& problem,
                                  IntersectionIterator& nit,
                                  const DomainType& x,
                                  RangePairType& pEn,
                                  RangePairType& pNeigh,
                                  ResultType& resEn, 
                                  ResultType& resNeigh) {
      problem.numericalFlux(nit, x, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            Element<2>::get(pEn), Element<2>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType,
      class IntersectionIterator,
      class QuadratureType,
      class ResultType
      >
    static void numericalFlux(ProblemType& problem,
                              IntersectionIterator& nit,
                              QuadratureType& quad,
                              int quadPoint,
                              RangePairType& pEn,
                              RangePairType& pNeigh,
                              ResultType& resEn,
                              ResultType& resNeigh) {
      problem.numericalFlux(nit, quad, quadPoint, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            Element<2>::get(pEn), Element<2>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       const DomainType& x,
                       RangePairType& p,
                       ResultType& res) {
      problem.source(en, x, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     Element<2>::get(p), 
                     res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void source(ProblemType& problem,
                       QuadratureType& quad,
                       Entity& en,
                       int quadPoint,
                       RangePairType& p,
                       ResultType res) {
      problem.source(en, quad, quadPoint, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     Element<2>::get(p), 
                     res);
    }
  };

  // Four arguments
  template <class T1, class T2, class T3, class T4>
  class Caller<Pair<T1, Pair<T2, Pair<T3, Pair<T4, Nil> > > > > {
  public:
    typedef Pair<T1, Pair<T2, Pair<T3, Pair<T4, Nil> > > > RangePairType;
    
    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               const DomainType& x,
                               RangePairType& p, 
                               ResultType& res) {
      problem.analyticalFlux(en, x, 
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             Element<2>::get(p), 
                             Element<3>::get(p), 
                             res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               QuadratureType& quad,
                               int quadPoint,
                               RangePairType& p,
                               ResultType& res) {
      problem.analyticalFlux(en, quad, quadPoint,
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             Element<2>::get(p), 
                             Element<3>::get(p), 
                             res);
    }

    template <
      class ProblemType, 
      class IntersectionIterator, 
      class DomainType, 
      class ResultType
      >
    static void numericalFluxFlux(ProblemType& problem,
                                  IntersectionIterator& nit,
                                  const DomainType& x,
                                  RangePairType& pEn,
                                  RangePairType& pNeigh,
                                  ResultType& resEn, 
                                  ResultType& resNeigh) {
      problem.numericalFlux(nit, x, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            Element<2>::get(pEn), Element<2>::get(pNeigh),
                            Element<3>::get(pEn), Element<3>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType,
      class IntersectionIterator,
      class QuadratureType,
      class ResultType
      >
    static void numericalFlux(ProblemType& problem,
                              IntersectionIterator& nit,
                              QuadratureType& quad,
                              int quadPoint,
                              RangePairType& pEn,
                              RangePairType& pNeigh,
                              ResultType& resEn,
                              ResultType& resNeigh) {
      problem.numericalFlux(nit, quad, quadPoint, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            Element<2>::get(pEn), Element<2>::get(pNeigh),
                            Element<3>::get(pEn), Element<3>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       const DomainType& x,
                       RangePairType& p,
                       ResultType& res) {
      problem.source(en, x, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     Element<2>::get(p), 
                     Element<3>::get(p), 
                     res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       QuadratureType& quad,
                       int quadPoint,
                       RangePairType& p,
                       ResultType res) {
      problem.source(en, quad, quadPoint, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     Element<2>::get(p), 
                     Element<3>::get(p), 
                     res);
    }
  };

  // Five argument
  template <class T1, class T2, class T3, class T4, class T5>
  class Caller<Pair<T1, Pair<T2, Pair<T3, Pair<T4, Pair<T5, Nil> > > > > > {
  public:
    typedef Pair<T1, Pair<T2, Pair<T3, Pair<T4, Pair<T5, Nil> > > > > RangePairType;
    
    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               const DomainType& x,
                               RangePairType& p, 
                               ResultType& res) {
      problem.analyticalFlux(en, x, 
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             Element<2>::get(p), 
                             Element<3>::get(p), 
                             Element<4>::get(p), 
                             res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void analyticalFlux(ProblemType& problem,
                               Entity& en,
                               QuadratureType& quad,
                               int quadPoint,
                               RangePairType& p,
                               ResultType& res) {
      problem.analyticalFlux(en, quad, quadPoint,
                             Element<0>::get(p), 
                             Element<1>::get(p), 
                             Element<2>::get(p), 
                             Element<3>::get(p), 
                             Element<4>::get(p), 
                             res);
    }

    template <
      class ProblemType, 
      class IntersectionIterator, 
      class DomainType, 
      class ResultType
      >
    static void numericalFluxFlux(ProblemType& problem,
                                  IntersectionIterator& nit,
                                  const DomainType& x,
                                  RangePairType& pEn,
                                  RangePairType& pNeigh,
                                  ResultType& resEn, 
                                  ResultType& resNeigh) {
      problem.numericalFlux(nit, x, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            Element<2>::get(pEn), Element<2>::get(pNeigh),
                            Element<3>::get(pEn), Element<3>::get(pNeigh),
                            Element<4>::get(pEn), Element<4>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType,
      class IntersectionIterator,
      class QuadratureType,
      class ResultType
      >
    static void numericalFlux(ProblemType& problem,
                              IntersectionIterator& nit,
                              QuadratureType& quad,
                              int quadPoint,
                              RangePairType& pEn,
                              RangePairType& pNeigh,
                              ResultType& resEn,
                              ResultType& resNeigh) {
      problem.numericalFlux(nit, quad, quadPoint, 
                            Element<0>::get(pEn), Element<0>::get(pNeigh),
                            Element<1>::get(pEn), Element<1>::get(pNeigh),
                            Element<2>::get(pEn), Element<2>::get(pNeigh),
                            Element<3>::get(pEn), Element<3>::get(pNeigh),
                            Element<4>::get(pEn), Element<4>::get(pNeigh),
                            resEn, resNeigh);
    }

    template <
      class ProblemType, class Entity, class DomainType, class ResultType>
    static void source(ProblemType& problem,
                       Entity& en,
                       const DomainType& x,
                       RangePairType& p,
                       ResultType& res) {
      problem.source(en, x, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     Element<2>::get(p), 
                     Element<3>::get(p), 
                     Element<4>::get(p), 
                     res);
    }

    template <
      class ProblemType, class Entity, class QuadratureType, class ResultType>
    static void source(ProblemType& problem,
                       QuadratureType& quad,
                       Entity& en,
                       int quadPoint,
                       RangePairType& p,
                       ResultType res) {
      problem.source(en, quad, quadPoint, 
                     Element<0>::get(p), 
                     Element<1>::get(p), 
                     Element<2>::get(p), 
                     Element<3>::get(p), 
                     Element<4>::get(p), 
                     res);
    }
  };
*/
  
} // end namespace Dune
  
#endif
