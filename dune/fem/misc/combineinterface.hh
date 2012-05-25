#ifndef DUNE_COMBINEINTERFACE_HH
#define DUNE_COMBINEINTERFACE_HH

#include <dune/fem/misc/femtuples.hh>

namespace Dune {
/*! @addtogroup InterfaceHelper
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

    @include example.cc

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

    At the time being we only allow combinations of reference types
    to interfaces since copying is difficult for BN
 */
template <class T1,class T2>  
struct PairOfInterfaces {
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1 t1,T2 t2); 
  //! Access first interface instance
  T1Type& first();
  //! Access second interface instance
  T2Type& second();
  //! const Access first interface instance
  const T1Type& first() const ;
  //! const Access second interface instance
  const T2Type& second() const;
};
template <class T1,class T2>  
struct PairOfInterfaces<T1&,T2&> {  
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1& t1,T2& t2) : p2(t2,n), p(t1,p2) {}
  //! Access first interface instance
  T1Type& first() {return p.first();}
  //! Access second interface instance
  T2Type& second() {return p.second().first();}
  //! const Access first interface instance
  const T1Type& first() const {return p.first();}
  //! const Access second interface instance
  const T2Type& second() const {return p.second().first();}
private:
  typedef Pair<T2&,Nil> Pair2Type;
  typedef Pair<T1&,Pair2Type> PairType;
  Nil n;
  Pair2Type p2;
  PairType p;
};
template <class T1,class T2>  
struct PairOfInterfaces<T1*,T2&> {  
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1* t1,T2& t2) : p2(t2,n), p(t1,p2) {}
  //! Access first interface instance
  T1Type& first() {return *(p.first());}
  //! Access second interface instance
  T2Type& second() {return p.second().first();}
  //! Access first interface instance
  const T1Type& first() const {return *(p.first());}
  //! Access second interface instance
  const T2Type& second() const {return p.second().first();}
private:
  typedef Pair<T2&,Nil> Pair2Type;
  typedef Pair<T1*,Pair2Type> PairType;
  Nil n;
  Pair2Type p2;
  PairType p;
};
template <class T1,class T2>  
struct PairOfInterfaces<T1&,T2*> {  
  //! First interface classes as non-reference types
  typedef typename TypeTraits<T1>::ReferredType T1Type;
  //! Second interface classes as non-reference types
  typedef typename TypeTraits<T2>::ReferredType T2Type;
  //! Constructor taking two interface instances
  PairOfInterfaces(T1& t1,T2* t2) : p2(t2,n), p(t1,p2) {}
  //! Access first interface instance
  T1Type& first() {return p.first();}
  //! Access second interface instance
  T2Type& second() {return *(p.second().first());}
  //! Access first interface instance
  const T1Type& first() const {return p.first();}
  //! Access second interface instance
  const T2Type& second() const {return *(p.second().first());}
private:
  typedef Pair<T2*,Nil> Pair2Type;
  typedef Pair<T1&,Pair2Type> PairType;
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
  T1Type& first() {return *(p.first());}
  //! Access second interface instance
  T2Type& second() {return *(p.second().first());}
  //! Access first interface instance
  const T1Type& first() const {return *(p.first());}
  //! Access second interface instance
  const T2Type& second() const {return *(p.second().first());}
private:
  typedef Pair<T2*,Nil> Pair2Type;
  typedef Pair<T1*,Pair2Type> PairType;
  Nil n;
  Pair2Type p2;
  PairType p;
};
/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
// Helper class for combining Operators
template<template <class T11,class T21> class CI,
	 typename T1, typename T2, typename T3, typename T4, typename T5,
	 typename T6, typename T7, typename T8, typename T9>
struct TupleToInterfacePair {
  typedef CI<T1,typename 
	     TupleToInterfacePair<CI,T2,T3,T4,T5,T6,T7,T8,T9,Nil>::Type&> Type;
};
template<template <class T11,class T21> class CI,
	 typename T1,typename T2>
struct TupleToInterfacePair<CI,T1,T2,Nil,Nil,Nil,Nil,Nil,Nil,Nil> {
  typedef CI<T1,T2> Type;
};
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
struct CombineInterface : 
  public TupleToInterfacePair<CI,T1,T2,T3,T4,T5,T6,T7,T8,T9>::Type {
  typedef CombineInterface<CI,T1,T2,T3,T4,T5,T6,T7,T8,T9> ThisType;
  typedef CombineInterface<CI,T2,T3,T4,T5,T6,T7,T8,T9,Nil> NextType;
  typedef typename TupleToInterfacePair<CI,T1,T2,T3,T4,T5,T6,T7,T8,T9>::Type
  BaseType;
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2,
		   typename TupleAccessTraits<T3>::ParameterType t3,
		   typename TupleAccessTraits<T4>::ParameterType t4,
		   typename TupleAccessTraits<T5>::ParameterType t5,
		   typename TupleAccessTraits<T6>::ParameterType t6,
		   typename TupleAccessTraits<T7>::ParameterType t7,
		   typename TupleAccessTraits<T8>::ParameterType t8,
		   typename TupleAccessTraits<T9>::ParameterType t9)
    : BaseType(t1,next) 
    , next(t2,t3,t4,t5,t6,t7,t8,t9) {}
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2,
		   typename TupleAccessTraits<T3>::ParameterType t3,
		   typename TupleAccessTraits<T4>::ParameterType t4,
		   typename TupleAccessTraits<T5>::ParameterType t5,
		   typename TupleAccessTraits<T6>::ParameterType t6,
		   typename TupleAccessTraits<T7>::ParameterType t7,
		   typename TupleAccessTraits<T8>::ParameterType t8)
    : BaseType(t1,next) 
    , next(t2,t3,t4,t5,t6,t7,t8) {}
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2,
		   typename TupleAccessTraits<T3>::ParameterType t3,
		   typename TupleAccessTraits<T4>::ParameterType t4,
		   typename TupleAccessTraits<T5>::ParameterType t5,
		   typename TupleAccessTraits<T6>::ParameterType t6,
		   typename TupleAccessTraits<T7>::ParameterType t7)
    : BaseType(t1,next) 
    , next(t2,t3,t4,t5,t6,t7) {}
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2,
		   typename TupleAccessTraits<T3>::ParameterType t3,
		   typename TupleAccessTraits<T4>::ParameterType t4,
		   typename TupleAccessTraits<T5>::ParameterType t5,
		   typename TupleAccessTraits<T6>::ParameterType t6)
    : BaseType(t1,next) 
    , next(t2,t3,t4,t5,t6) {}
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2,
		   typename TupleAccessTraits<T3>::ParameterType t3,
		   typename TupleAccessTraits<T4>::ParameterType t4,
		   typename TupleAccessTraits<T5>::ParameterType t5)
    : BaseType(t1,next) 
    , next(t2,t3,t4,t5) {}
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2,
		   typename TupleAccessTraits<T3>::ParameterType t3,
		   typename TupleAccessTraits<T4>::ParameterType t4)
    : BaseType(t1,next) 
    , next(t2,t3,t4) {}
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2,
		   typename TupleAccessTraits<T3>::ParameterType t3)
    : BaseType(t1,next) 
    , next(t2,t3) {}
  CombineInterface(const ThisType& ci)
    : BaseType(ci.first(),next) 
    , next(ci.next) {}
  NextType next;
};
template <template <class T11,class T21> class CI,class T1,class T2> 
struct CombineInterface<CI,T1,T2,Nil,Nil,Nil,Nil,Nil,Nil,Nil>
  : public CI<T1,T2> {
 public:
  typedef CI<T1,T2> BaseType;
  typedef CombineInterface<CI,T1,T2,Nil,Nil,Nil,Nil,Nil,Nil,Nil> 
  ThisType;  
  CombineInterface(typename TupleAccessTraits<T1>::ParameterType t1,
		   typename TupleAccessTraits<T2>::ParameterType t2)
    : BaseType(t1,t2) {}
  CombineInterface(const ThisType& ci) 
    : BaseType(ci.first(),ci.second()) {}
};
// *******************************************************
// *******************************************************

} // end namespace Dune 
#endif
