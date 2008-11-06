#ifndef DUNE_SELECTION_HH
#define DUNE_SELECTION_HH

#include <dune/common/misc.hh>
#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/misc/utility.hh> 
#include <dune/fem/pass/pass.hh>

namespace Dune
{

  /**
   * @brief Converts a selector definition to a set of pairs.
   *
   * A selector is a tuple where the integer constants are converted using
   * Int2Type and the end marker is mapped from -1 to Nil.
   */
  template <int N1, int N2, int N3, int N4, int N5, int N6, int N7, int N8, int N9>
  struct SelectorToPairs {
    typedef Pair<
      Int2Type<N1>, 
      typename SelectorToPairs<N2, N3, N4, N5, N6, N7, N8, N9, -1>::Type
    > Type;
  };

  /*
   * @brief Specialisation for the closure.
   *
   * For selectors, -1 serves as the termination marker.
   */
  template <>
  struct SelectorToPairs<-1, -1, -1, -1, -1, -1, -1, -1, -1> {
    typedef Nil Type;
  };

  /**
   * @brief A list of compile-time constants to select entries of a tuple.
   */

  template <int N1 = -1,
            int N2 = -1, 
            int N3 = -1, 
            int N4 = -1, 
            int N5 = -1, 
            int N6 = -1, 
            int N7 = -1, 
            int N8 = -1, 
            int N9 = -1>
  struct Selector : 
    public SelectorToPairs<N1, N2, N3, N4, N5, N6, N7, N8, N9> 
  {
    typedef typename SelectorToPairs<
      N1, N2, N3, N4, N5, N6, N7, N8, N9>::Type Base;
  private:
    //! A Selector contains pure type information, so there is no need to have
    //! instances
    Selector();
  };
  

  /**
   * @brief Gives back passNum of the pass which passId
   * is equal to template-given id. Also checks if all 
   * passes have passId
   *
   */
  template< class Pass, int id >
  struct FindPass
  {
    CompileTimeChecker< (Pass::passId != -1) > __ASSERT_ALL_PASSES_HAVE_PASSID_OR_NONE__;
    typedef typename Pass :: PreviousPassType PreviousPassType;
    enum { passNum = ( (int)Pass::passId == (int)id ? (int)Pass :: passNum
                     : (int)FindPass< PreviousPassType, id > :: passNum ) };
  private:
    // no need to have instance of this class
    FindPass();
  };


  template < class Argument , int StartPassIdImp, int id >
  struct FindPass< StartPass< Argument , StartPassIdImp > , id >
  {
    enum { passNum = ( (int)StartPass< Argument , StartPassIdImp >
                        ::passId == (int)id ? 0 : -1 ) };
  };


  /**
   *
   * @brief Converts passId to passDiff ( number which tells how much to count
   * backwards from the last pass to get the pass which passId is
   * equal to template-given id, with exception that 0 means StartPass )
   *
   */
  template < class Pass , int id >
  struct PassId2PassDiff
  {
    enum { passNum = (int)FindPass< typename Pass::PreviousPassType , id >::passNum };
    CompileTimeChecker< (passNum != -1) > __ASSERT_PASS_ID_FOUND_AND_ASSERT_ALL_PASSES_HAVE_ID_OR_NONE;
    enum { passDiff = ((passNum == 0) ? 0 : Pass::passNum - passNum) };
  };

  template< class Pass , int id , bool passHasId >
  struct CompatiblePassId2PassDiff;

  template< class Pass , int id >
  struct CompatiblePassId2PassDiff< Pass , id , true >
  { 
    // in this case template-given id is passId
    enum { passDiff = PassId2PassDiff< Pass , id >::passDiff };
  };
  
  template< class Pass , int id >
  struct CompatiblePassId2PassDiff< Pass , id , false >
  { 
    // in this case template-given id is already passDiff
    enum { passDiff = id };
    // check if all pases don't have passId
    CompileTimeChecker< Pass::PreviousPassType::passId == -1 > __ASSERT_ALL_PASSES_HAVE_PASSID_OR_NONE__;
    typedef CompatiblePassId2PassDiff< typename Pass::PreviousPassType , id , false > CheckPassIds;
  };
  

  /**
   *
   * @brief Carries PassType and Selector to Filter
   *
   */
  template< class Pass , class Selector > 
  struct CombinedSelector;
  

  template< class Pass , int N1, int N2, int N3, int N4, int N5, int N6, int N7, int N8, int N9 > 
  struct CombinedSelector< Pass 
    , Selector< N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9 > >
  {
    typedef Pass PassType;
    typedef typename Selector< N1 , N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9 >::Base Base;
    typedef CombinedSelector< Pass , Selector< N2 , N3 , N4 , N5 , N6 , N7 , N8 , N9 , -1 > > Type2;
  };
  
  

  template <class SelectorType>
  struct MaxIndex {
    enum { value = -1 };
  };

  template <class Head, class Tail>
  struct MaxIndex<Pair<Head, Tail> > {
    enum { value = (static_cast<int>(Head::value) > 
                    static_cast<int>(MaxIndex<Tail>::value) ?
                    static_cast<int>(Head::value) : 
                    static_cast<int>(MaxIndex<Tail>::value)) };
  };

  template <>
  struct MaxIndex<Nil> {
    enum { value = -1 };
  };

  
  /**
   *
   * @brief For checking if everything is ok with types
   *
   */
  template< class > struct SelectorPrint;

  template< class H , class T>
  struct SelectorPrint< Pair< H , T > >
  {
    static void apply()
    {
      std::cout <<H::value <<"\n";
      SelectorPrint< T >::apply();
    }
  };

  template< >
  struct SelectorPrint< Nil >
  {
    static void apply() {}
  };

  
  /**
   *
   * @brief Converts TupleId to TupleDiff ( number which tells how much to count
   * backwards from the last tuple element to get the tuple element which value is
   * equal to template-given Id )
   *
   */
  template< class SelectorBaseImp , int id , int diff = 0 >
  struct Id2convertedId
  {
    enum { num = ( (int)SelectorBaseImp::Type1::value == id ? diff :
                   Id2convertedId< typename SelectorBaseImp::Type2 , id , diff+1 >::num ) };
    //CompileTimeChecker< (num != -1) > __ASSERT_ARGUMENT_TUPLE_ID_FOUND__;
  };

  // specialization 
  template< int id , int diff >
  struct Id2convertedId< Selector< -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 , -1 >::Base , id , diff >
  {
    enum { num =  -1 };
  };


  /**
   *
   * @brief MatchTuplesForType gives the type of Tuple
   * entry which has same offset in tuple like the Selector entry 
   * that matches the template-given id
   *
   */
  template< class Pair , class Tuple , int id , bool matched >
  struct MatchTuplesHelperForType
  {};
  
  template< class SelectorHead , class SelectorTail , 
            class TupleHead , class TupleTail , int id >
  struct MatchTuplesHelperForType< Pair< SelectorHead , SelectorTail > , 
                            Pair< TupleHead , TupleTail > , id , false >
  {
    typedef typename MatchTuplesHelperForType
      < SelectorTail , TupleTail, id 
        , (int)SelectorTail::Type1::value == id >::Type Type;
  };
  
  template< class SelectorHead ,  
            class TupleHead , class TupleTail , int id >
  struct MatchTuplesHelperForType< Pair< SelectorHead , Nil > , 
                                   Pair< TupleHead , TupleTail > , id , false >
  {
    // when unable to match selector values with template-given id
    // then return Nil as a type
    typedef Nil Type;
  };

  template< class SelectorHead , class SelectorTail , 
            class TupleHead , class TupleTail , int id >
  struct MatchTuplesHelperForType< Pair< SelectorHead , SelectorTail > , 
                            Pair< TupleHead , TupleTail > , id , true >
  {
    typedef TupleHead Type;
  };

  template< class SelectorBase , class Tuple , int id >
  struct MatchTuplesForType
  {
    typedef typename MatchTuplesHelperForType
      < SelectorBase , Tuple , id , (int)SelectorBase::Type1::value == id >
      ::Type Type;
  };



  /**
   *
   * @brief MatchTuplesForValue gives the value of Tuple
   * entry which has same offset in tuple like the Selector entry 
   * that matches the template-given id
   *
   */
  template< class Pair, class Tuple, int id, bool match >
  struct MatchTuplesHelperForValue
  {};
  
  template< class SelectorHead, class SelectorTail, 
            class TupleHead, class TupleTail, int id >
  struct MatchTuplesHelperForValue
    < Pair< SelectorHead, SelectorTail >, Pair< TupleHead, TupleTail >, id, false >
  {
    typedef Pair< SelectorHead, SelectorTail > Selector;
    typedef Pair< TupleHead, TupleTail > Tuple;

  private:
    typedef typename MatchTuplesForType< Selector, Tuple, id > :: Type MatchedType;

    static const bool matchTail = ((int)SelectorTail :: Type1 :: value == id);

  public:
    static typename TupleAccessTraits< MatchedType > :: ConstType
    get( const Tuple &arg )
    {
      return MatchTuplesHelperForValue< SelectorTail, TupleTail, id, matchTail >
        :: get( arg.second() );
    }
    
    static typename TupleAccessTraits< MatchedType > :: NonConstType
    get( Tuple &arg )
    {
      return MatchTuplesHelperForValue< SelectorTail , TupleTail, id, matchTail >
        :: get( arg.second() );
    }
  };
  
  template< class SelectorHead, class TupleHead, class TupleTail, int id >
  struct MatchTuplesHelperForValue
    < Pair< SelectorHead, Nil >, Pair< TupleHead, TupleTail >, id, false >
  {
    typedef Pair< SelectorHead, Nil > Selector;
    typedef Pair< TupleHead, TupleTail > Tuple;

    // when unable to match selector values with template-given id
    // then return nulltype() as a value
    static typename TupleAccessTraits< Nil > :: ConstType
    get( const Tuple &arg )
    {
      return Nil();
    }
    
    static typename TupleAccessTraits< Nil > :: NonConstType
    get( Tuple &arg )
    {
      return Nil();
    }
  };

  template< class SelectorHead, class SelectorTail,
            class TupleHead, class TupleTail, int id >
  struct MatchTuplesHelperForValue
    < Pair< SelectorHead, SelectorTail >, Pair< TupleHead, TupleTail >, id, true >
  {
    typedef Pair< SelectorHead, SelectorTail > Selector;
    typedef Pair< TupleHead, TupleTail > Tuple;

    static typename TupleAccessTraits< TupleHead > :: ConstType
    get( const Tuple &arg )
    {
      return arg.first();
    }
    
    static typename TupleAccessTraits< TupleHead > :: NonConstType
    get( Tuple &arg )
    {
      return arg.first();
    }
  };

  template< class SelectorBase, class Tuple, int id >
  struct MatchTuplesForValue
  {
  private:
    typedef typename MatchTuplesForType< SelectorBase, Tuple, id > :: Type MatchedType;

    static const bool match = ((int)SelectorBase::Type1::value == id);

  public:
    static typename TupleAccessTraits< MatchedType > :: ConstType
    get( const Tuple &arg )
    {
      return MatchTuplesHelperForValue< SelectorBase, Tuple, id, match > :: get( arg );
    }
    
    static typename TupleAccessTraits< MatchedType > :: NonConstType
    get( Tuple &arg )
    {
      return MatchTuplesHelperForValue< SelectorBase, Tuple, id, match > :: get( arg );
    }
  };



  template< class Selector, class Head, class Tail >
  struct SelectorPair
  : public Pair< Head, Tail >
  {
    typedef SelectorPair< Selector, Head, Tail > ThisType;
    typedef Pair< Head, Tail > BaseType;

  public:
    // For compatibility with ElementAccess, ElementType
    typedef BaseType FirstPair;

    typedef Selector SelectorType;
    
    template< int id >
    struct Get
    {
      //enum { convertedId = Id2convertedId< SelectorType , id >::num };
      //typedef typename ElementType< convertedId, BaseType > :: Type Type;
      //typedef typename ElementType< id, BaseType > :: Type Type;
      typedef typename MatchTuplesForType< SelectorType , BaseType , id > :: Type Type;
    };

  public:
    inline SelectorPair ( typename TupleAccessTraits< Head > :: ParameterType head,
                          Tail &tail )
    : BaseType( head, tail )
    {}

    inline SelectorPair ( const ThisType &other )
    : BaseType( other )
    {}

    inline explicit SelectorPair ( const BaseType &other )
    : BaseType( other )
    {}
    
    template< int id >
    inline const typename Get< id > :: Type &get () const
    {
      //enum { convertedId = Id2convertedId< SelectorType , id >::num };
      //return Element< convertedId > :: get( (const BaseType&)(*this) );
      //return Element< id > :: get( (const BaseType&)(*this) );
      return MatchTuplesForValue< SelectorType , BaseType , id > :: get( (const BaseType&)(*this) );
    }

    template< int id >
    inline typename Get< id > :: Type &get ()
    {
      //enum { convertedId = Id2convertedId< SelectorType , id >::num };
      //return Element< convertedId > :: get( (BaseType&)(*this) );
      //return Element< id > :: get( (BaseType&)(*this) );
      return MatchTuplesForValue< SelectorType , BaseType , id > :: get( (BaseType&)(*this) );
    }

    template< int id >
    inline const typename Get< id > :: Type &
    operator[] ( const Int2Type< id > idVariable ) const
    {
      return get< id >();
    }

    template< int id >
    inline typename Get< id > :: Type &
    operator[] ( const Int2Type< id > idVariable )
    {
      return get< id >();
    }
  };

} // end namespace Dune

#endif
