#ifndef DUNE_PHONYIIT
#define DUNE_PHONYIIT

namespace Dune
{
  
  template< class IntersectionType, class IntersectionIteratorType >
  struct PhonyIntersectionIterator
  {
    typedef IntersectionIteratorType Type;
  };
  
  template< class IntersectionIteratorType >
  struct PhonyIntersectionIterator< IntersectionIteratorType, IntersectionIteratorType >
  {
    typedef int Type;
  };
 
}

#endif
