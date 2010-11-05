#ifndef DUNE_CALLERUTILITY_HH
#define DUNE_CALLERUTILITY_HH

#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/misc/utility.hh>
#include <dune/common/typetraits.hh>

#include <dune/fem/function/localfunction/temporarylocalfunction.hh>

#include "selection.hh"

namespace Dune
{

  /**
   * @brief Extracts elements from a tuple based on a Selector definition
   *                                                                                                                                                                           
   * Use this class to extract the elements from a tuple. The selector 
   * specifies the elements to be extracted by indicating their index. The
   * extracted elements are again stored in a tuple.
   */
  template< class ArgTupleImp, class SelectorImp >
  struct Filter 
  {
    dune_static_assert( ((int)MaxIndex< SelectorImp >::value < (int)TupleLength< ArgTupleImp >::value),
                        "Maximum index of selector exceeds argument length." );
    
    typedef typename SelectorImp::Base SelectorType;
    typedef typename SelectorImp::PassType PassType;

    static const bool hasPassId = ((int)SelectorImp::PassType::passId != -1);

    //! The index of the extracted element
    static const int index
      = CompatiblePassId2PassDiff< PassType, ElementType< 0, SelectorType >::Type::value, hasPassId >::passDiff;

    //! The type of the extracted element.
    typedef typename ElementType< index, ArgTupleImp >::Type AppendType;

    //! The type of the next stage of the filter.
    typedef Filter< ArgTupleImp, typename SelectorImp::Type2 > NextFilterType;

    typedef typename NextFilterType::ResultType NextResultType;
    
    //! The filtered tuple resulting from this filter stage.
    // typedef Pair<AppendType, typename NextFilterType::ResultType> ResultType;
    typedef SelectorPair< SelectorType, AppendType, NextResultType > ResultType;
    
    //! Extracts the elements specified by the selector template argument.
    //! \param arg Argument tuple.
    static ResultType apply ( ArgTupleImp &arg )
    {
      typedef Pair< AppendType, NextResultType > BasicResultType;
      NextResultType nextResult = NextFilterType::apply( arg );
      return ResultType( BasicResultType( Element< index >::get( arg ), nextResult ) );
    }
  };
  
  /**
   * @brief Specialisation for closure.
   *
   * The end of the filtering is marked by the last element of the selector.
   * \warning If the largest index in the selector is greater than the length
   *          of the tuple, an ugly compile-time error will occur.
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
    // old version uses local function of discrete function
    //typedef typename TypeTraits<DFType>::PointeeType::LocalFunctionType Type;
    typedef typename TypeTraits<DFType>::PointeeType DiscreteFunctionType;
    typedef ConstLocalFunction< DiscreteFunctionType > Type;
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
    typedef typename TypeTraits<Head>::PointeeType DiscreteFunctionType;
    typedef ConstLocalFunction< DiscreteFunctionType > LocalFunctionType;

  public:
    static inline ResultType apply(const Pair<Head, Tail>& pairs) {
      LocalFunctionType tmp( *pairs.first() );
      //LocalFunctionType tmp(pairs.first()->localFunctionStorage());
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
    //typedef typename TypeTraits<Head>::PointeeType::LocalFunctionType LocalFunctionType;
    typedef typename TypeTraits<Head>::PointeeType DiscreteFunctionType ;
    //typedef typename DiscreteFunctionType :: LocalFunctionType  LocalFunctionType;
    // old version uses local function of discrete function
    typedef ConstLocalFunction< DiscreteFunctionType > LocalFunctionType;
    typedef Pair<LocalFunctionType, Nil> ResultType;
 
  public:
    static inline ResultType apply(const Pair<Head, Nil>& pairs) {
      //LocalFunctionType tmp(pairs.first()->localFunctionStorage());
      LocalFunctionType tmp(*pairs.first());
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
    LocalFunctionSetter( const EntityImp& en ) :
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
    const EntityImp& en_;
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

  template <class VectorTupleType, int passId >
  class TupleToVectorConverter
  {

    TupleToVectorConverter(const TupleToVectorConverter&);
  public:
    typedef typename VectorTupleType :: value_type TupleType;
    typedef typename ElementType<passId, TupleType> :: Type ValueType;

    TupleToVectorConverter(VectorTupleType& vec)
      : vector_( vec )
    {}

    ValueType& operator [] (const size_t i)
    {
      assert( i < vector_.size() );
      return Element< passId > :: get( vector_[ i ] );
    }

    const ValueType& operator [] (const size_t i) const
    {
      assert( i < vector_.size() );
      return Element< passId > :: get( vector_[ i ] );
    }

    size_t size() const 
    {
      return vector_.size();
    }

  protected:  
    VectorTupleType& vector_;
  };

  template <class TupleType1, class VectorType>
  class ForEachValueVector {
  public:
    //! Constructor
    //! \param t1 First tuple.
    //! \param t2 Second tuple.
    ForEachValueVector(TupleType1& t1, VectorType& vec ) :
      tuple1_(t1),
      vector_( vec )
    {}

    //! Applies the function object f to the pair of tuples.
    //! \param f The function object to apply on the pair of tuples.
    template <class Functor>
    void apply(Functor& f) 
    {
      apply <0> (f, tuple1_);
    }

  private:
    //! Specialisation for the last element.
    template <int passId, class Functor, class Head1>
    void apply(Functor& f, Pair<Head1, Nil>& last1) 
    {
      callFunction< passId > ( f, last1.first() ); 
    }

    //! Specialisation for a standard element.
    template <int passId, class Functor, class Head1, class Tail1>
    void apply(Functor& f, Pair<Head1, Tail1>& p1) 
    {
      callFunction< passId > ( f, p1.first() ); 
      apply< passId + 1 >( f, p1.second() );
    }

    template <int passId, class Functor, class P1> 
    void callFunction(Functor& f, P1& p1) 
    {
      TupleToVectorConverter< VectorType, passId > vector ( vector_ ); 
      f.visit( p1, vector );
    }

  private:
    TupleType1& tuple1_;
    VectorType& vector_;
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
  class LocalFunctionEvaluateQuadrature 
  {
  public:
    //! Constructor
    //! \param quad The quadrature in question.
    LocalFunctionEvaluateQuadrature(const QuadratureImp& quad )
      : quad_( quad )
    {}

    //! Evaluation of a local function
    template <class LFType, class RangeVectorType>
    void visit(LFType& lf, RangeVectorType& ranges) 
    {
      lf.evaluateQuadrature ( quad_ , ranges);
    }

  private:
    LocalFunctionEvaluateQuadrature();
    LocalFunctionEvaluateQuadrature(const LocalFunctionEvaluateQuadrature&);
    LocalFunctionEvaluateQuadrature& operator=(const LocalFunctionEvaluateQuadrature&);

  private:
    const QuadratureImp& quad_;
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
