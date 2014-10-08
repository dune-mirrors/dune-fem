#ifndef DUNE_FEM_COMMON_COORDINATE_HH
#define DUNE_FEM_COMMON_COORDINATE_HH

namespace Dune
{

  namespace Fem
  {

    template< class Point >
    static inline const Point &coordinate ( const Point &x )
    {
      return x;
    }

  } // namespace Fem

} // namespace Dune

#endif // #ifndef DUNE_FEM_COMMON_COORDINATE_HH
