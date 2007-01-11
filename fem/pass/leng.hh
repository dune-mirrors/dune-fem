#define DUNE_LEN_HH
#ifndef DUNE_LEN_HH

#include<ostream>

namespace Dune 
{
#  ifdef BOOST_NO_INCLASS_MEMBER_INITIALIZATION
#       define BOOST_STATIC_CONSTANT(type, assignment) enum { assignment }
#  else
#     define BOOST_STATIC_CONSTANT(type, assignment) static const type assignment
#  endif

template<class T>
struct length  {
  BOOST_STATIC_CONSTANT(int, value = 1 + length<typename T::Type2>::value);
};

template<>
 struct length<Tuple<Nil,Nil,Nil,Nil,Nil,Nil,Nil,Nil,Nil> > {
  BOOST_STATIC_CONSTANT(int, value = 0);
};

template<>
struct length<Nil> {
  BOOST_STATIC_CONSTANT(int, value = 0);
};
  
}

#endif
