#include <dune/common/tuples.hh>
namespace Dune {
/*! @defgroup CombineInterface Interface helpers
    A general strategy for combining up to nine instances of an 
    interface class. The resulting type again satisfies the 
    interface so it can again be combined for further instances or
    combined objects. 
    The first template argument must be a template class which
    combines two instances of an interface and itself satisfies the 
    interface - this can be compared with the Pair class as basis of
    the Tuple class. Derive this class from the interface and the
    helper class PairOfInterfaces<T1,T2>.
    The other nine template arguments should be inmplemetations
    of the interface. 
    
    Note that to avoid copying of instances reference or pointer
    types should be used. 

    @example @include /home/andreas/dune/testing/example.cc

    Further usage showing for example how to use types defined
    in the interface class and how to automatically
    combine these to new interface classes can be found in...
 
 */
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3=Nil,class T4=Nil,
	  class T5=Nil,class T6=Nil,class T7=Nil,class T8=Nil,class T9=Nil> 
struct CombineInterface;
/*! Use as base class for all InterfacePair classes. Handles
    storage of interface objects and defines
    first() and second() method to access the two interface objects.
 */
template <class T1,class T2>  
struct PairOfInterfaces {  
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1 t1,T2 t2) : p2(t2,n), p(t1,p2) {}
  //! Access first interface instance
  T1Type& first() {return p.first();}
  //! Access second interface instance
  T2Type& second() {return p.second().first();}
private:
  typedef Pair<T2,Nil> Pair2Type;
  typedef Pair<T1,Pair2Type> PairType;
  Nil n;
  Pair2Type p2;
  PairType p;
};
template <class T1,class T2>  
struct PairOfInterfaces<T1*,T2> {  
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1* t1,T2 t2) : p2(t2,n), p(t1,p2) {}
  //! Access first interface instance
  T1Type& first() {return *(p.first());}
  //! Access second interface instance
  T2Type second() {return p.second().first();}
private:
  typedef Pair<T2,Nil> Pair2Type;
  typedef Pair<T1*,Pair2Type> PairType;
  Nil n;
  Pair2Type p2;
  PairType p;
};
template <class T1,class T2>  
struct PairOfInterfaces<T1,T2*> {  
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1 t1,T2* t2) : p2(t2,n), p(t1,p2) {}
  //! Access first interface instance
  T1Type& first() {return p.first();}
  //! Access second interface instance
  T2Type second() {return *(p.second().first());}
private:
  typedef Pair<T2*,Nil> Pair2Type;
  typedef Pair<T1,Pair2Type> PairType;
  Nil n;
  Pair2Type p2;
  PairType p;
};
template <class T1,class T2>  
struct PairOfInterfaces<T1*,T2*> {  
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1* t1,T2* t2) : p2(t2,n), p(t1,p2) {}
  //! Access first interface instance
  T1Type first() {return *(p.first());}
  //! Access second interface instance
  T2Type second() {return *(p.second().first());}
private:
  typedef Pair<T2*,Nil> Pair2Type;
  typedef Pair<T1*,Pair2Type> PairType;
  Nil n;
  Pair2Type p2;
  PairType p;
};
// Helper class for combining Operators
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4,
	  class T5,class T6,class T7,class T8,class T9>
struct CombMultInterHelper;
/*!  Wrapper class to combine interface objects to a single interface object.

     Template arguments: the class for combining two single interfaces and
                         a list of interface classes which should be combined
                         (from two to nine).

     The resulting class is derived from the combined interface class CI
     and is therefore itself an interface implementation.
*/
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4,
	  class T5,class T6,class T7,class T8,class T9>
struct CombineInterface : public CombMultInterHelper<CI,T1,T2,T3,T4,T5,T6,T7,T8,T9> {
  CombineInterface(T1 t1,T2 t2,T3 t3=Nil(),T4 t4=Nil(),
		   T5 t5=Nil(),T6 t6=Nil(),T7 t7=Nil(),T8 t8=Nil(),T9 t9=Nil()) 
    : CombMultInterHelper<CI,T1,T2,T3,T4,T5,T6,T7,T8,T9>(t1,t2,t3,t4,t5,t6,t7,t8,t9)
  {}
};
/********************************************************************/
// The helper classes which take care of different number of
// template parameters
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4,
	  class T5,class T6,class T7,class T8,class T9>
struct CombMultInterHelper : 
    public CI<T1,CombMultInterHelper<CI,T2,T3,T4,T5,T6,T7,T8,T9,Nil> >  {
  typedef CombMultInterHelper<CI,T2,T3,T4,T5,T6,T7,T8,T9,Nil> NextType;
  typedef CI<T1,NextType> BaseType;
  CombMultInterHelper(T1 t1,T2 t2,T3 t3,T4 t4,T5 t5,T6 t6,T7 t7,T8 t8,T9 t9) :
    BaseType(t1,NextType(t2,t3,t4,t5,t6,t7,t8,t9,Nil())) {}
}; 
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4,
	  class T5,class T6,class T7,class T8>
struct CombMultInterHelper<CI,T1,T2,T3,T4,T5,T6,T7,T8,Nil> : 
    public CI<T1,CombMultInterHelper<CI,T2,T3,T4,T5,T6,T7,T8,Nil,Nil> >  {
  typedef CombMultInterHelper<CI,T2,T3,T4,T5,T6,T7,T8,Nil,Nil> NextType;
  typedef CI<T1,NextType> BaseType;
  CombMultInterHelper(T1 t1,T2 t2,T3 t3,T4 t4,T5 t5,T6 t6,T7 t7,T8 t8,
		      Nil n9) :
    BaseType(t1,NextType(t2,t3,t4,t5,t6,t7,t8,Nil(),Nil())) {}
}; 
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4,
	  class T5,class T6,class T7>
struct CombMultInterHelper<CI,T1,T2,T3,T4,T5,T6,T7,Nil,Nil> : 
    public CI<T1,CombMultInterHelper<CI,T2,T3,T4,T5,T6,T7,Nil,Nil,Nil> >  {
  typedef CombMultInterHelper<CI,T2,T3,T4,T5,T6,T7,Nil,Nil,Nil> NextType;
  typedef CI<T1,NextType> BaseType;
  CombMultInterHelper(T1 t1,T2 t2,T3 t3,T4 t4,T5 t5,T6 t6,T7 t7,
		      Nil n8,Nil n9) :
    BaseType(t1,NextType(t2,t3,t4,t5,t6,t7,Nil(),Nil(),Nil())) {}
}; 
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4,
	  class T5,class T6>
struct CombMultInterHelper<CI,T1,T2,T3,T4,T5,T6,Nil,Nil,Nil> : 
    public CI<T1,CombMultInterHelper<CI,T2,T3,T4,T5,T6,Nil,Nil,Nil,Nil> >  {
  typedef CombMultInterHelper<CI,T2,T3,T4,T5,T6,Nil,Nil,Nil,Nil> NextType;
  typedef CI<T1,NextType> BaseType;
  CombMultInterHelper(T1 t1,T2 t2,T3 t3,T4 t4,T5 t5,T6 t6,
		      Nil n7,Nil n8,Nil n9) :
    BaseType(t1,NextType(t2,t3,t4,t5,t6,Nil(),Nil(),Nil(),Nil())) {}
}; 
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4,
	  class T5>
struct CombMultInterHelper<CI,T1,T2,T3,T4,T5,Nil,Nil,Nil,Nil> : 
    public CI<T1,CombMultInterHelper<CI,T2,T3,T4,T5,Nil,Nil,Nil,Nil,Nil> >  {
  typedef CombMultInterHelper<CI,T2,T3,T4,T5,Nil,Nil,Nil,Nil,Nil> NextType;
  typedef CI<T1,NextType> BaseType;
  CombMultInterHelper(T1 t1,T2 t2,T3 t3,T4 t4,T5 t5,
		      Nil n6,Nil n7,Nil n8,Nil n9) :
    BaseType(t1,NextType(t2,t3,t4,t5,Nil(),Nil(),Nil(),Nil(),Nil())) {}
}; 
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3,class T4>
struct CombMultInterHelper<CI,T1,T2,T3,T4,Nil,Nil,Nil,Nil,Nil> : 
    public CI<T1,CombMultInterHelper<CI,T2,T3,T4,Nil,Nil,Nil,Nil,Nil,Nil> >  {
  typedef CombMultInterHelper<CI,T2,T3,T4,Nil,Nil,Nil,Nil,Nil,Nil> NextType;
  typedef CI<T1,NextType> BaseType;
  CombMultInterHelper(T1 t1,T2 t2,T3 t3,T4 t4,
		      Nil n5,Nil n6,Nil n7,Nil n8,Nil n9) :
    BaseType(t1,NextType(t2,t3,t4,Nil(),Nil(),Nil(),Nil(),Nil(),Nil())) {}
}; 
template <template <class T11,class T21> class CI,
	  class T1,class T2,class T3>
struct CombMultInterHelper<CI,T1,T2,T3,Nil,Nil,Nil,Nil,Nil,Nil> : 
    public CI<T1,CombMultInterHelper<CI,T2,T3,Nil,Nil,Nil,Nil,Nil,Nil,Nil> >  {
  typedef CombMultInterHelper<CI,T2,T3,Nil,Nil,Nil,Nil,Nil,Nil,Nil> NextType;  
  typedef CI<T1,NextType> BaseType;
  CombMultInterHelper(T1 t1,T2 t2,T3 t3,
		      Nil n4,Nil n5,Nil n6,Nil n7,Nil n8,Nil n9) :
    BaseType(t1,NextType(t2,t3,Nil(),Nil(),Nil(),Nil(),Nil(),Nil(),Nil())) {}
}; 
template <template <class T11,class T21> class CI,
	  class T1,class T2>
struct CombMultInterHelper<CI,T1,T2,Nil,Nil,Nil,Nil,Nil,Nil,Nil> : public CI<T1,T2>  {
  CombMultInterHelper(T1 t1,T2 t2,
		      Nil n3,Nil n4,Nil n5,Nil n6,Nil n7,Nil n8,Nil n9) :  
    CI<T1,T2>(t1,t2) {}
};
// *******************************************************
// *******************************************************
}
