#ifndef DUNE_CALLERUTILITY_HH
#define DUNE_CALLERUTILITY_HH

#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/common/typetraits.hh>

#include "selection.hh"

namespace Dune {

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
    typedef CompileTimeChecker<
      static_cast<int>(MaxIndex<SelectorImp>::value) < 
      static_cast<int>(TupleLength<ArgTupleImp>::value)
    > Maximal_index_of_selector_exceeds_argument_length;
    
    //typedef SelectorImp SelectorType;
    typedef typename SelectorImp::Base SelectorType;

    //! The index of the extracted element.
    enum { index = (int)CompatiblePassId2PassDiff< typename SelectorImp::PassType
                     , ElementType<0, SelectorType>::Type::value 
                     , ((int)SelectorImp::PassType::passId != -1) > :: passDiff };

    //! The type of the extracted element.
    typedef typename ElementType<index, ArgTupleImp>::Type AppendType;

    //! The type of the next stage of the filter.
    typedef Filter< ArgTupleImp , typename SelectorImp::Type2 > NextFilterType;
    
    //! The filtered tuple resulting from this filter stage.
    // typedef Pair<AppendType, typename NextFilterType::ResultType> ResultType;
    typedef SelectorPair
      < SelectorType, AppendType, typename NextFilterType :: ResultType >
      ResultType;
    
    //! Extracts the elements specified by the selector template argument.
    //! \param arg Argument tuple.
    static inline ResultType apply(ArgTupleImp& arg) {
      typename NextFilterType::ResultType nextFilter = NextFilterType::apply(arg);
      return ResultType( Pair< AppendType , typename NextFilterType::ResultType >
                         ( Element<index>::get(arg) , nextFilter ) );
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

  template < class ArgTupleImp , class Pass , int N1>
  struct Filter<ArgTupleImp, CombinedSelector< Pass , Selector<N1,-1,-1,-1,-1,-1,-1,-1,-1> > >
  {
    typedef typename Selector<N1,-1,-1,-1,-1,-1,-1,-1,-1>::Base SelectorType;

    enum { index = (int)CompatiblePassId2PassDiff< Pass , N1 
                   , (Pass::passId != -1) > :: passDiff };
    typedef typename ElementType< index , ArgTupleImp>::Type AppendType;

    //! The closure of a tuple is always Nil
    typedef SelectorPair< SelectorType , AppendType , Nil > ResultType;
    
    //! Provides the closure of the filtered tuple
    static inline ResultType apply(ArgTupleImp& arg) {
      return ResultType( Pair< AppendType , Nil >(Element<index>::get(arg) , nullType() ) );
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
   * @brief Extracts the type of the Range vector from a tuple of 
   * LocalFunction pointers.
   */
  template <class LFType>
  struct RangeTypeEvaluator {
    typedef typename LFType::RangeType Type;
  };

  /**
   * @brief Extracts the type of the jacobian range from a tuple of 
   * LocalFunction pointers.
   */
  template <class LFType>
  struct JacobianRangeTypeEvaluator {
    typedef typename LFType::JacobianRangeType Type;
  };


  // * Idea: can I do the same with ForEachValue?
  // * Give tuple in constructor and get the rest with apply()
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
      return Result(LocalFunction(*p.first()),nullType());
    }

    template <class Head, class Tail>
    static inline typename Traits<Head, Tail>::ResultType
    apply(Pair<Head, Tail>& p) {
      typedef typename Traits<Head, Tail>::ResultType Result;
      typedef typename Traits<Head, Tail>::LocalFunctionType LocalFunction;
      return Result(LocalFunction(*p.first()), 
                    LocalFunctionCreator2::apply(p.second()));
    }
  };
  
  /**
   * @brief Creates a tuple of local function pointers from a tuple of 
   * discrete function pointers by calling newLocalFunction on the discrete
   * function (plus a copy constructor...)
   * \note Don't forget to delete the local functions in the end!
   */
  template< class DFTupleType >
  class LocalFunctionCreator
  {};

  template< class Selector, class Head, class Tail >
  class LocalFunctionCreator< SelectorPair< Selector, Head, Tail > >
  : public LocalFunctionCreator< Pair< Head, Tail > >
  {
    typedef LocalFunctionCreator< Pair< Head, Tail > > BaseType;

  public:
    typedef typename BaseType :: ResultType PairResultType;
    typedef SelectorPair
      < Selector, typename PairResultType :: Type1, typename PairResultType :: Type2 >
      ResultType;

    static inline ResultType apply( const Pair< Head, Tail > &pairs )
    {
      return ResultType( BaseType :: apply( pairs ) );
    }
  };

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
    static inline ResultType apply(const Pair<Head, Tail>& pairs) {
      LocalFunctionType tmp(pairs.first()->localFunctionStorage());
      typename LocalFunctionCreator<Tail>::ResultType localFunc = LocalFunctionCreator<Tail>::apply(pairs.second());
      return ResultType(tmp,localFunc);
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
    static inline ResultType apply(const Pair<Head, Nil>& pairs) {
      LocalFunctionType tmp(pairs.first()->localFunctionStorage());
      return ResultType(tmp, nullType());
    }
  };

  /**
   * @brief Generates a tuple based on types defined by another tuple.
   */
  template <template <class> class TypeEvaluator, class LFTupleType>
  class Creator
  {};

  template< template< class > class TypeEvaluator,
            class Selector, class Head, class Tail >
  class Creator< TypeEvaluator, SelectorPair< Selector, Head, Tail > >
  : public Creator< TypeEvaluator, Pair< Head, Tail > >
  {
    typedef Creator< TypeEvaluator, Pair< Head, Tail > > BaseType;

  public:
    typedef typename BaseType :: ResultType PairResultType;
    typedef SelectorPair
      < Selector, typename PairResultType :: Type1, typename PairResultType :: Type2 >
      ResultType;

    static inline ResultType apply ()
    {
      return ResultType( BaseType :: apply() );
    }
    
    static inline ResultType apply ( Pair< Head, Tail > &pairs )
    {
      return ResultType( BaseType :: apply( pairs ) );
    }
  };
  
  /**
   * @brief Specialisation for standard element.
   */
  template <template <class> class TypeEvaluator, class Head, class Tail>
  class Creator<TypeEvaluator, Pair<Head, Tail> > {
  private:
    template <template <class> class TE, class T>
    friend class Creator;
    typedef typename ForEachType<TypeEvaluator, Tail>::Type TailType;
  
  public:
    typedef typename ForEachType<
      TypeEvaluator, Pair<Head, Tail> >::Type ResultType;
    typedef typename TypeEvaluator<Head>::Type ValueType;

  public:
    static inline ResultType apply() {
      ValueType tmp(0.0);
      typename Creator<TypeEvaluator, Tail>::ResultType created = Creator<TypeEvaluator, Tail>::apply();
      ResultType result(tmp, created);
      return result;
    }

    static inline ResultType apply(Pair<Head, Tail>& pairs) {
      ValueType tmp(0.0);
      typename Creator<TypeEvaluator, Tail>::ResultType created = Creator<TypeEvaluator, Tail>::apply(pairs.second());
      return ResultType(tmp,created);
    }
  };
  
  /**
   * @brief Specialisation for last element
   */
  template <template <class> class TypeEvaluator, class Head>
  class Creator<TypeEvaluator, Pair<Head, Nil> > {
    template <template <class> class TE, class T>
    friend class Creator;

  public:
    typedef typename TypeEvaluator<Head>::Type ValueType;
    typedef Pair<ValueType, Nil> ResultType;

  public:
    static inline ResultType apply() {
      ValueType tmp(0.0);
      return ResultType(tmp, nullType());
    }

    static inline ResultType apply(Pair<Head, Nil>& pairs) {
      ValueType tmp(0.0);
      return ResultType(tmp, nullType());
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
    template <class LFType>
    void visit(LFType& lf) {
       lf.init(en_);
    }

  private:
    LocalFunctionSetter();
    LocalFunctionSetter(const LocalFunctionSetter&);
    LocalFunctionSetter& operator=(const LocalFunctionSetter&);

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
  template <class DomainImp>
  class LocalFunctionEvaluateLocal {
  public:
    //! Constructor
    //! \param x The local coordinate on en.
    LocalFunctionEvaluateLocal(const DomainImp& x) :
      x_(x)
    {}

    //! Triggers the evaluation of a local function
    template <class LFType, class RangeType>
    void visit(LFType& lf, RangeType& res) {
      lf.evaluate(x_, res);
    }

  private:
    LocalFunctionEvaluateLocal();
    LocalFunctionEvaluateLocal(const LocalFunctionEvaluateLocal&);
    LocalFunctionEvaluateLocal& operator=(const LocalFunctionEvaluateLocal&);
    
  private:
    const DomainImp& x_;
  };

  /**
   * @brief  Calls lf.evaluate on a tuple of local functions with a 
   * corresponding tuple of results.
   *
   * Use this Functor in conjunction with a ForEachValuePair built from a tuple
   * of local functions and a tuple with the corresponding range vectors.
   *
   */
  template <class QuadratureImp>
  class LocalFunctionEvaluateQuad {
    typedef typename QuadratureImp :: QuadraturePointWrapperType QuadraturePointWrapperType;
  public:
    //! Constructor
    //! \param quad The quadrature in question.
    //! \param quadPoint The index of the quadrature point of quadrature quad
    LocalFunctionEvaluateQuad(const QuadratureImp& quad,
                              const int quadPoint) :
      quadPoint_( quad[quadPoint] )
    {}

    //! Evaluation of a local function
    template <class LFType, class RangeType>
    void visit(LFType& lf, RangeType& res) {
      lf.evaluate( quadPoint_ , res);
    }

  private:
    LocalFunctionEvaluateQuad();
    LocalFunctionEvaluateQuad(const LocalFunctionEvaluateQuad&);
    LocalFunctionEvaluateQuad& operator=(const LocalFunctionEvaluateQuad&);

  private:
    const QuadraturePointWrapperType quadPoint_;
  };

  template <class DomainImp>
  class LocalFunctionEvaluateJacobianLocal {
  public:
    LocalFunctionEvaluateJacobianLocal(const DomainImp& x) :
      x_(x)
    {}

    template <class LFType, class JRangeType>
    void visit(LFType& lf, JRangeType& res) {
      lf.jacobian(x_, res);
    }

  private:
    LocalFunctionEvaluateJacobianLocal();
    LocalFunctionEvaluateJacobianLocal(const LocalFunctionEvaluateJacobianLocal&);
    LocalFunctionEvaluateJacobianLocal& 
    operator=(const LocalFunctionEvaluateJacobianLocal&);

  private:
    const DomainImp& x_;
  };
  
  template <class QuadratureImp>
  class LocalFunctionEvaluateJacobianQuad {
    typedef typename QuadratureImp :: QuadraturePointWrapperType QuadraturePointWrapperType;
  public:
    LocalFunctionEvaluateJacobianQuad(const QuadratureImp& quad,
                                      const int quadPoint) :
      quadPoint_ ( quad[quadPoint] )
    //  quad_(quad),
    //  quadPoint_(quadPoint)
    {}

    template <class LFType, class JRangeType>
    void visit(LFType& lf, JRangeType& res) {
      //lf.jacobian(quad_, quadPoint_, res);
      lf.jacobian(quadPoint_, res);
    }

  private:
    LocalFunctionEvaluateJacobianQuad();
    LocalFunctionEvaluateJacobianQuad(const LocalFunctionEvaluateJacobianQuad&);
    LocalFunctionEvaluateJacobianQuad& 
    operator=(const LocalFunctionEvaluateJacobianQuad&);

  private:
    const QuadraturePointWrapperType quadPoint_;
    //const QuadratureImp& quad_;
    //const int quadPoint_;
  };

  /**
   * @brief  Calls lf.setQuad on a tuple of local function.
   *
   * Use this Functor in conjunction with a ForEachValue built from a tuple
   * of local functions.
   *
   */
  template <class EntityImp, class QuadratureImp>
  class DiscreteFunctionSetQuad {
  public:
    //! Constructor
    //! \param en the entity
    //! \param quad The quadrature in question.
    DiscreteFunctionSetQuad(EntityImp& en,QuadratureImp& quad) :
      quad_(quad),
      en_(en)
    {}

    //! Set the quadrature for a local function
    template <class LFType>
    void visit(LFType& lf) {
    }

  private:
    DiscreteFunctionSetQuad();
    DiscreteFunctionSetQuad(const DiscreteFunctionSetQuad&);
    DiscreteFunctionSetQuad& operator=(const DiscreteFunctionSetQuad&);

  private:
    QuadratureImp& quad_;
    EntityImp& en_;
  };

  template <class QuadratureImp>
  class LocalDiscreteFunctionSetQuad {
  public:
    //! Constructor
    //! \param quad The quadrature in question.
    LocalDiscreteFunctionSetQuad(QuadratureImp& quad) :
      quad_(quad)
    {}

    //! Set the quadrature for a local function
    template <class LFType>
    void visit(LFType& lf) {
     lf.getBaseFunctionSet().addQuadrature(quad_);
    }

  private:
    LocalDiscreteFunctionSetQuad();
    LocalDiscreteFunctionSetQuad(const LocalDiscreteFunctionSetQuad&);
    LocalDiscreteFunctionSetQuad& operator=(const LocalDiscreteFunctionSetQuad&);

  private:
    QuadratureImp& quad_;
  };

} // end namespace Dune
#endif
