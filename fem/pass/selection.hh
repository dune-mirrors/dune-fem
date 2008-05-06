#ifndef DUNE_SELECTION_HH
#define DUNE_SELECTION_HH

#include <dune/fem/misc/femtuples.hh>
#include <dune/fem/misc/utility.hh> 
#include <dune/common/misc.hh>

namespace Dune {
  /**
   * @brief Converts a selector definition to a set of pairs.
   *
   * A selector is a tuple where the integer constants are converted using
   * Int2Type and the end marker is mapped from -1 to Nil.
   */
  template <int N1, 
            int N2, 
            int N3, 
            int N4, 
            int N5, 
            int N6, 
            int N7, 
            int N8, 
            int N9>
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
  

  // go through pass tree and give pass number
  // of which passId is equal to template-given id
  template< class Pass, int id >
  struct FindPass
  {
    CompileTimeChecker< (Pass::passId != -1) > __USE_ONLY_ONE_TYPE_OF_SELECTOR_IN_ALL_DMODELS__;
    //typedef typename Pass :: DiscreteModelType DMType;
    typedef typename Pass :: PreviousPassType PreviousPassType;
    //enum {  passNum = ( (int)DMType::id == (int)id ? (int)Pass :: passNum
    //                 : (int)FindPass< PreviousPassType, id > :: passNum ) };
    enum {  passNum = ( (int)Pass::passId == (int)id ? (int)Pass :: passNum
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
     * @brief Given a pass tree and a selector, that takes ids of passes
     * it gives a selector that contains corresponding pass numbers.
     *
     */
  template< class Pass, class Selector >
  struct ConvertSelector;


  // switch from type Selector to type SelectorToPoints = Selector :: Base
  template< class Pass, int N1, int N2, int N3, int N4, int N5, int N6, int N7, int N8, int N9 >
  struct ConvertSelector< Pass, Selector< N1, N2, N3, N4, N5, N6, N7, N8, N9 > >
  : public ConvertSelector< Pass, typename Selector< N1, N2, N3, N4, N5, N6, N7, N8, N9 > :: Base >
  {};


  template< class Pass, int id, class Tail >
  struct ConvertSelector< Pass, Pair< Int2Type< id >, Tail > >
  {
    enum { passNum = FindPass< typename Pass::PreviousPassType , id > :: passNum };
    enum { passDiff = ((passNum == 0) ? 0 : Pass::passNum - passNum) };
    typedef Pair< Int2Type< passDiff >, typename ConvertSelector< Pass, Tail > :: Base > Base;

    // in case pass with model.id = id hasn't been founded
    // then produce compile error ( by creating instance of undefined
    // class CompileTimeChecker< false >
    CompileTimeChecker< (passNum != -1) > __Assert_Pass_Id_Found__;
  private:
    // no need to have instance of this class
    ConvertSelector();
  };


  template< class Pass >
  struct ConvertSelector< Pass, Nil >
  {
    typedef Nil Base;
  };


#if 1
  template< class Pass , class Selector , bool old >
  struct CompatibleConvertSelector;

  
  template< class Pass , class Selector >
  struct CompatibleConvertSelector< Pass , Selector , true >
  {
    typedef typename Selector::Base Base;
  private:
    // no need to have instance of this class
    CompatibleConvertSelector();
  };

  
  template< class Pass , class Selector >
  struct CompatibleConvertSelector< Pass , Selector , false >
  {
    typedef typename ConvertSelector< Pass , Selector >::Base Base;
  };
#endif


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

} // end namespace Dune

#endif
